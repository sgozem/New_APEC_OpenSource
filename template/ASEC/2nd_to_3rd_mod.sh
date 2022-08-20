#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
tempdir=`grep "tempdir" ../Infos.dat | awk '{ print $2 }'`

opt=${Project}_VDZ_Opt

if [[ -f $opt/$opt.out ]]; then
   echo ""
   echo " *************************************************************"
   echo "                      Warning!"
   echo ""
   echo " $opt.out already excist. We are goint to use it..."
   echo ""
   echo " *************************************************************"
   echo ""
else
   echo ""
   echo " Collecting the CASSCF/VDZ optimization from iRODS..."
   echo ""
   dir=`basename $tempdir`
   iget -r /arctic/projects/CHEM9C4/$USER/$dir $opt
   if [[ -f $opt/$dir/$opt.out ]]; then
      mv $opt/$dir/* $opt
      rm -r $opt/$dir
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

#
# Instructions to the user
#
echo ""
echo " The current project is $Project. Checking the CASSCF/VDZ optimization..."
echo ""

#
# Checking if the optimization ended successfully, or if the assigned hours ended
#
opt=${Project}_VDZ_Opt
if [ -d $opt ]; then
   cd $opt
   if grep -q "Timing: Wall=" $opt.out; then
      echo " CAS optimization ended successfully"
      echo ""
   else
      echo " The CAS optimization did not finished well."
      echo " Please check if the assigned hours expired or if some error occurred."
      echo ""
      echo " If you need to restart, run the restart.sh script."
      echo " Otherwise find out what is the error and re-run 1st_to_2nd.sh."
      echo ""
      echo " 2nd_to_3rd.sh will terminate now!"
      echo ""
      cp $templatedir/restart.sh ../
      exit 0
   fi
else
   echo " CAS/VDZ optimization folder was not found! Check what is wrong"
   echo " Terminating"
   exit 0 
fi

#
# If the folder exists, the script is aborted with an error
#
cd ..
new=${Project}_VDZP
if [ -d $new ]; then
   ./smooth_restart.sh $new "Do you want to re-run the QM/MM VDZP single point? (y/n)" 7
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir ${new}
cp $opt/$opt.key ${new}/${new}.key
cp $opt/$opt.Final.xyz ${new}/${new}.xyz
cp $opt/${prm}.prm ${new}/
cp $templatedir/molcas.slurm.sh ${new}/molcas-job.sh
cp $templatedir/ASEC/templateSP ${new}/
cd ${new}/

#
# Editing the template for the CAS single point
#
sed -i "s|PARAMETRI|${prm}|" templateSP

mv templateSP ${new}.input
sed -i "s/ANO-L-VDZ/ANO-L-VDZP/g" ${new}.input

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
sed -i "s|hh:00:00|60:00:00|" molcas-job.sh

#
# Job submission and template copy for the following step
#
echo ""
echo " Submitting the CAS/VDZP single point now..."
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
cp $templatedir/ASEC/sp_to_opt_VDZP_mod.sh .
cp $templatedir/ASEC/alter_orbital_mod.sh .
../update_infos.sh "Next_script" "sp_to_opt_VDZP_mod.sh" ../Infos.dat
echo ""
echo "******************************************************************"
echo ""
echo " As soon as the calculation is finished, run sp_to_opt_VDZP_mod.sh"
echo ""
echo "******************************************************************"
echo ""

