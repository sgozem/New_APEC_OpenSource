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
echo " The current project is $Project. Checking the CAS/6-31G* optimization..."
echo ""

# Checking if the optimization ended successfully, or if the assigned hours ended
#
opt=${Project}_6-31G_Opt
if [ -d $opt ]; then
   cd $opt
   contr=`grep 'landing' $opt.out | awk '{ print $1 }'`
   if [[ $contr == "Happy" ]]; then
      echo " CAS/6-31G* optimization ended successfully"
      echo ""
   else
      echo " The CAS/6-31G* optimization did not finished well."
      echo " Please check if the assigned hours expired or if some error occurred."
      echo ""
      echo " If you need to restart, run restart.sh"
      echo " Otherwise find out what is the error and re-run sp_to_opt_631G.sh."
      echo ""
      echo " 3rd_to_4th.sh will terminate now!"
      echo ""
      cp $templatedir/restart.sh ../
      exit 0
   fi
else
   echo " CAS/6-31G* optimization folder was not found! Check what is wrong"
   echo " Terminating"
   exit 0 
fi

# Retrieving the occupation numbers from CAS single point for checking purposes
#
grep -A2 "Natural orbitals and occupation numbers for root  1" ${opt}.out

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
            echo " Going ahead with the CASPT2/6-31G* single point"
	    echo ""
	 else
	    echo " You might have a problem with active space selection. To fix it:"
            echo " 1. - Use Molden to look at the orbitals in $opt.rasscf.molden"
	    echo " 2. - Find the right orbital(s) to be placed in the active space"
	    echo " 3. - Run the script alter_orbital.sh"
	    echo ""
	    echo " Now 3rd_to_4th.sh will terminate"
	    echo ""
	    cp $templatedir/alter_orbital.sh .
	    exit 0
	 fi
      fi
done

# Preparing the files for CASPT2/6-31G* single point
# If the folder exists, the script is aborted with an error
#
cd ..
new=${Project}_CASPT2
if [ -d $new ]; then
   ./smooth_restart.sh $new "Do you want to re-run the CASPT2 single point? (y/n)" 9
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir ${new}
cp $opt/$opt.JobIph ${new}/${new}.JobIph
cp $opt/$opt.key ${new}/${new}.key
cp $opt/$opt.Final.xyz ${new}/${new}.xyz
cp $opt/${prm}.prm ${new}/
cp $templatedir/molcas-job.sh ${new}/
cp $templatedir/templateCaspt2 ${new}/
cp $templatedir/modify-inp.vim ${new}/
cd ${new}/

# Editing the template for the CASPT2 single point
#
sed -i "s|PARAMETRI|${prm}|" templateCaspt2

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

# Merging the geometrical part and the template for CASPT2 single point
#
mv ${new}.input temp
cat temp templateCaspt2 > ${new}.input
rm temp templateCaspt2

# Writing project name, input directory and memory requested in the submission script
# A CASPT2 needs 3500 MB and 24 hrs...
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
sed -i "s|MEMTOT|3800|" molcas-job.sh
sed -i "s|MEMORIA|3500|" molcas-job.sh
sed -i "s|hh:00:00|24:00:00|" molcas-job.sh

# Job submission and template copy for the following step
#
echo ""
echo " Submitting the CASPT2 single point now..."
echo ""
sleep 1

qsub molcas-job.sh

cd ..
cp $templatedir/finalPDB.sh .
cp $templatedir/resultcaspt2.sh .

# Updating the Infos.dat Current field, which stores the current running calculation
#
awk '{ if ( $1 == "CurrCalc") sub($2,"CASPT2"); print}' ../Infos.dat > temp
mv temp ../Infos.dat

# Messages to the user
#
echo ""
echo ""
echo " As soon as the calculation is finished, you can use:"
echo " - finalPDB.sh, to generate the PDB file of the structure"
echo " - resultcaspt2.sh, to calculate excitation energies and absorption maxima"
echo ""

