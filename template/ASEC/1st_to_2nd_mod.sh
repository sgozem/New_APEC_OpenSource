#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
tempdir=`grep "tempdir" ../Infos.dat | awk '{ print $2 }'`

sp=${Project}_VDZ

if [[ -f $sp/$sp.out ]]; then
   echo ""
   echo " *************************************************************"
   echo "                      Warning!"
   echo ""
   echo " $sp.out already excist. We are goint to use it..."
   echo ""
   echo " *************************************************************"
   echo ""
else
   echo ""
   echo " Collecting the SP CASSCF/VDZ from iRODS..."
   echo ""
   dir=`basename $tempdir`
   iget -r /arctic/projects/CHEM9C4/$USER/$dir $sp
   if [[ -f $sp/$dir/$sp.out ]]; then
      mv $sp/$dir/* $sp
      rm -r $sp/$dir
      irm -r /arctic/projects/CHEM9C4/$USER/$dir
   else
      echo ""
      echo "************************************************************************"
      echo ""
      echo " It seems that the MD is still running or it did not finish properly"
      echo ""
      echo "************************************************************************"
      echo ""
      exit 0
   fi
fi

# Instructions to the user
#
echo ""
echo " The current project is $Project. Checking the SP CASSCF/VDZ single point..."
echo ""

#
# Checking if the single point ended successfully, with control on folder existence
#
if [ -d $sp ]; then
   cd $sp
   if grep -q "Timing: Wall=" $sp.out; then
      echo " CAS single point ended successfully" 
      echo ""
   else
      echo " CAS single point still in progress. Terminating..."
      echo ""
      exit 0
   fi
else
   echo " CAS/3-21G single point folder not found! Check what is wrong"
   echo " Terminating..."
   exit 0
fi

#
# Retrieving the occupation numbers from CAS single point for checking purposes
#
grep -A2 "Natural orbitals and occupation numbers for root  1" ${sp}.out
cd ..

#
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
            echo " Going ahead with the CAS/3-21G optimization"
	    echo ""
	 else
	    echo " You might have a problem with active space selection. To fix it:"
            echo " 1. - Use Molden to look at the orbitals in $sp.rasscf.molden"
	    echo " 2. - Find the right orbital(s) to be placed in the active space"
	    echo " 3. - Run the script alter_orbital_mod.sh"
	    echo ""
	    echo " Now 1st_to_2nd.sh will terminate"
	    echo ""
	    cp $templatedir/ASEC/alter_orbital_mod.sh .
	    exit 0
	 fi
      fi
done

#
# Preparing the files for CASSCF/VDZ optimization 
# If the folder exists, the script is aborted with an error
#
new=${Project}_VDZ_Opt
if [ -d $new ]; then
   ./smooth_restart.sh $new "Do you want to re-run the QM/MM 3-21G optimization? (y/n)" 6
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
cp $templatedir/molcas.slurm.sh ${new}/molcas-job.sh
cp $templatedir/ASEC/template_CASSCF_min ${new}/template
#cp $templatedir/modify-inp.vim ${new}/
cd ${new}/

#
# Editing the template for the CAS optimization
#
sed -i "s|PARAMETRI|${prm}|" template

mv template ${new}.input

actspace=`grep "Active_space" ../../Infos.dat | awk '{ print $2 }'`
if [[ $actspace == "custom" ]]; then
   actelec=`grep "Active_electrons" ../../Infos.dat | awk '{ print $2 }'`
   actorb=`grep "Active_orbitals" ../../Infos.dat | awk '{ print $2 }'`
   inactorb=`grep "Inactive_orbitals" ../../Infos.dat | awk '{ print $2 }'`

   sed -i "s|nActEl=10 0 0|nActEl=$actelec 0 0|" $new.input
   sed -i "s|Ras2=10|Ras2=$actorb|" $new.input
   sed -i "s|Inactive=62|Inactive=$inactorb|" $new.input
fi

sed -i "s|NOMEPROGETTO|${new}|" molcas-job.sh
no=$PWD
sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|23000|" molcas-job.sh
sed -i "s|MEMORIA|20000|" molcas-job.sh
sed -i "s|hh:00:00|120:00:00|" molcas-job.sh

#
# Job submission and template copy for the following step
#
echo ""
echo " Submitting the CASSCF/VDZ optimization now..."
echo ""
#sleep 1

#TMPFILE=`mktemp -d /scratch/photon_XXXXXX`
#../../update_infos.sh "tempdir" $TMPFILE ../../Infos.dat
#sed -i "s|TEMPFOLDER|$TMPFILE|g" molcas-job.sh
#cp -r * $TMPFILE
#current=$PWD
#cd $TMPFILE
sbatch molcas-job.sh
#cd $current

cd ..
cp $templatedir/ASEC/2nd_to_3rd_mod.sh .
../update_infos.sh "Next_script" "2nd_to_3rd_mod.sh" ../Infos.dat
echo ""
echo "**********************************************************"
echo ""
echo " As soon as the calculation is finished, run 2nd_to_3rd.sh"
echo ""
echo "**********************************************************"
echo ""

