#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
tempdir=`grep "tempdir" ../Infos.dat | awk '{ print $2 }'`

scf=${Project}_OptSCF

if [[ -f $scf/$scf.out ]]; then
   echo ""
   echo " *************************************************************"
   echo "                      Warning!"
   echo ""
   echo " $scf.out already excist. We are goint to use it..."
   echo ""
   echo " *************************************************************"
   echo ""
else
   echo ""
   echo " Collecting the HF optimization from iRODS..."
   echo ""
   dir=`basename $tempdir`
   iget -r /arctic/projects/CHEM9C4/$USER/$dir $scf
   if [[ -f $scf/$dir/$scf.out ]]; then
      mv $scf/$dir/* $scf
      rm -r $scf/$dir
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
echo " The current project is $Project. Checking the HF optimization..."
echo ""

#
# Grepping the Happy landing to check that the calculation ended up properly
#
if grep -q "Timing: Wall=" $scf/$scf.out; then
   echo " HF optimization ended successfully"
   echo ""
else
   echo " HF optimization still in progress. Terminating..."
   echo ""
   exit 0 
fi	

#
# Creation of the folder for CAS/3-21G single point and copy of all the files
# If the folder already exists, it finishes with an error message
#
new=${Project}_VDZ
if [ -d $new ]; then
   ./smooth_restart.sh $new "Do you want to re-run the QM/MM 3-21G single point? (y/n)" 5
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir ${new}
cp $scf/$scf.Final.xyz ${new}/${new}.xyz
cp $scf/$scf.key ${new}/${new}.key
cp $scf/${prm}.prm ${new}/
#cp $templatedir/modify-inp.vim ${new}/
#slurm
cp $templatedir/molcas.slurm.sh ${new}/molcas-job.sh
#cp $templatedir/molcas-job.sh ${new}/
cp $templatedir/ASEC/templateSP ${new}/
cd ${new}/

#
# Editing the template for single point
#
sed -i "s|PARAMETRI|${prm}|" templateSP

#
# Editing the submission script template for a CAS single point 
#

mv templateSP $new.input

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
# Submitting the CAS/3-21G single point
#
echo ""
echo " Submitting the CAS/ANO-L-VDZ single point now..."
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
cp $templatedir/ASEC/1st_to_2nd_mod.sh .
cp $templatedir/ASEC/alter_orbital_mod.sh .
../update_infos.sh "Next_script" "1st_to_2nd_mod.sh" ../Infos.dat

echo ""
echo "***************************************************************"
echo ""
echo " As soon as the CAS single point is finished, run 1st_to_2nd.sh"
echo ""
echo "***************************************************************"
echo ""


