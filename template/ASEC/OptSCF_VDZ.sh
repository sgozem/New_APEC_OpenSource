#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`

# Instructions to the user
#
echo ""
echo " The current project is $Project. Checking the HF optimization..."
echo ""

# Grepping the Happy landing to check that the calculation ended up properly
#
scf=${Project}_OptSCF
if grep -q "Timing: Wall=" $scf/$scf.out; then
   echo ""
   echo " HF/ANO-L-MB optimization ended successfully"
   echo ""
else
   echo " HF optimization still in progress. Terminating..."
   echo ""
   exit 0 
fi	

# Creation of the folder ${Project}_OptSCF_VDZ for single point and copy of all the files
# If the folder already exists, it finishes with an error message
#
new=${Project}_OptSCF_VDZ
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
cp $scf/${scf}.input ${new}/${new}.input
#slurm
cp $templatedir/molcas.slurm.sh ${new}/molcas-job.sh
#cp $templatedir/molcas-job.sh ${new}/
cd ${new}/


# Editing the submission script template for a CAS single point 
#
 
sed -i "s|NOMEPROGETTO|${new}|" molcas-job.sh
no=$PWD
sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|23000|" molcas-job.sh
sed -i "s|MEMORIA|20000|" molcas-job.sh
sed -i "s|hh:00:00|260:00:00|" molcas-job.sh
sed -i "s|Basis = ANO-L-MB|Basis = ANO-L-VDZ|g" ${new}.input

# Submitting the CAS/3-21G single point
#
echo ""
echo " Submitting the CAS/3-21G single point now..."
echo ""
sleep 1

sbatch molcas-job.sh

cd ..
cp $templatedir/ASEC/Molcami2_mod.sh .
../update_infos.sh "Next_script" "Molcami2_mod.sh" ../Infos.dat
echo ""
echo "***********************************************************"
echo ""
echo " As soon as the OptSCF_VDZ is finished, run Molcami2_mod.sh"
echo ""
echo "***********************************************************"


