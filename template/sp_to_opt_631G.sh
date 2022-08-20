#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
version=`grep "Tk" ../Infos.dat | awk '{ print $3 }'`

# Instructions to the user
#
echo ""
echo " The current project is $Project. Checking the CAS/6-31G* single point..."
echo ""

# Checking if the single point ended successfully, with control on folder existence
#
sp=${Project}_6-31G
if [ -d $sp ]; then
   cd $sp
   contr=`grep 'landing' $sp.out | awk '{ print $1 }'`
   if [[ $contr == "Happy" ]]; then
      echo " CAS/6-31G* single point ended successfully" 
      echo ""
   else
      echo " CAS/6-31G* single point still in progress. Terminating..."
      echo ""
      exit 0
   fi
else
   echo " CAS/6-31G* single point folder not found! Check what is wrong"
   echo " Terminating..."
   echo ""
   exit 0
fi

# Retrieving the occupation numbers from CAS single point for checking purposes
#
grep -A2 "Natural orbitals and occupation numbers for root  1" ${sp}.out

# Asking the user if the occupation numbers are ok
#
echo ""
echo -n " Are all the occupation numbers in the range 0.02 to 1.98? (y/n)"
read answer
echo ""
contr=0
while [  $contr = 0 ]; do
      if [[ $answer != "y" && $answer != "n" ]]; then
         echo -n " Please answer y or n... -> "
         read answer
	 echo ""
      else
	 contr=1
         if [[ $answer == "y" ]]; then
            echo " Going ahead with the CAS/6-31G* optimization"
	    echo ""
	 else
	    echo " You might have a problem with active space selection. To fix it:"
            echo " 1. - Use Molden to look at the orbitals in $sp.rasscf.molden"
	    echo " 2. - Find the right orbital(s) to be placed in the active space"
	    echo " 3. - Run the script alter_orbital.sh"
	    echo ""
	    echo " sp_to_opt_631G.sh will terminate"
	    echo ""
	    cp $templatedir/alter_orbital.sh .
	    exit 0
	 fi
      fi
done

# Preparing the files for CAS/6-31G* optimization 
# If the folder exists, the script is aborted with an error
#
cd ..
new=${Project}_6-31G_Opt
if [ -d $new ]; then
  ./smooth_restart.sh $new "Do you want to re-run the QM/MM 6-31G* optimization? (y/n)" 8
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir ${new}
cp $sp/$sp.JobIph ${new}/${new}.JobIph
cp $sp/$sp.key ${new}/${new}.key
cp $sp/$sp.xyz ${new}/${new}.xyz
cp $sp/${prm}.prm ${new}/
cp $sp/$sp.Espf.Data ${new}/${new}.Espf.Data
cp $templatedir/molcas-job.sh ${new}/
cp $templatedir/template ${new}/
cp $templatedir/modify-inp.vim ${new}/
cd ${new}/

# Editing the template for the CAS optimization
# Replacing Set maxiter with MOLCAS_MAXITER for Molcas 8
# Replacing PARAMETRI with current prm filename template
# 
#
case $version in
     4.2 | 5.1)
     sed -i "s|Set   maxiter   100|Set   maxiter   100|" template
     ;;
     6.2 | 6.3)
     sed -i "s|Set   maxiter   100|export MOLCAS_MAXITER = 100|" template
esac

sed -i "s|PARAMETRI|${prm}|" template

# Calling xyzedit for input generation
#
case $version in
     4.2 | 5.1)
     echo "20" >> tinkerxyzedit.vim
     echo "1" >> tinkerxyzedit.vim
     ;;
     6.2 | 6.3)
     echo "21" >> tinkerxyzedit.vim
     echo "1" >> tinkerxyzedit.vim
esac

$tinkerdir/xyzedit ${new}.xyz < tinkerxyzedit.vim
rm tinkerxyzedit.vim

# Editing the input file with modify-inp.vim: removal of the Tinker standard part
# and introduction of the basis set
#
sed -i 's/BASIS/6-31G*/' modify-inp.vim
sed -i 's/BAS2/6-31G\\*/' modify-inp.vim
vim -es $new.input < modify-inp.vim
rm modify-inp.vim

# Merging the geometrical part and the template for CAS optimization
#
mv ${new}.input temp
cat temp template > ${new}.input
rm temp template

# Writing the project name, the input directory, time and memory requested in the submission script
# Here there is a CASSCF/6-31G* optimization, so 2 Gb should be enough...
# The maximum available hours, 144 hrs (6 days), are requested
#
case $version in
     4.2 | 5.1)
     sed -i "s|MOLCASMEM=MEMORIAMB|MOLCASMEM=MEMORIAMB|" molcas-job.sh
     ;;
     6.2 | 6.3)
     sed -i "s|MOLCASMEM=MEMORIAMB|MOLCAS_MEM=MEMORIAMB|" molcas-job.sh
esac


sed -i "s|NOMEPROGETTO|${new}|" molcas-job.sh
no=$PWD
sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|2000|" molcas-job.sh
sed -i "s|MEMORIA|1800|" molcas-job.sh
sed -i "s|hh:00:00|144:00:00|" molcas-job.sh

# Job submission and template copy for the following step
#
echo ""
echo " Submitting the CAS/6-31G* optimization now..."
echo ""
sleep 1

qsub molcas-job.sh

cd ..
cp $templatedir/3rd_to_4th.sh .

# Updating the Infos.dat Current field, which stores the current running calculation
#
awk '{ if ( $1 == "CurrCalc") sub($2,"6-31G_Opt"); print}' ../Infos.dat > temp
mv temp ../Infos.dat

echo ""
echo ""
echo " As soon as the calculation is finished run 3rd_to_4th.sh"
echo ""
