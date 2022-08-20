#!/bin/bash
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" ../Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" ../Infos.dat | awk '{ print $2 }'`

mkdir CASPT2_ipea_025
mkdir CASPT2_ipea_0
#mkdir CASPT2_ipea_0/CASSCF_3_States

Project_new=${Project}_VDZP_Opt

cp $templatedir/ASEC/template_CASPT2 CASPT2_ipea_025/${Project}_CASPT2_025.input
cp $templatedir/ASEC/template_CASPT2 CASPT2_ipea_0/${Project}_CASPT2_0.input
sed -i "s/ipea = 0.25/ipea = 0.0/g" CASPT2_ipea_0/${Project}_CASPT2_0.input
sed -i "s|PARAMETER|${prm}|" CASPT2_ipea_025/${Project}_CASPT2_025.input
sed -i "s|PARAMETER|${prm}|" CASPT2_ipea_0/${Project}_CASPT2_0.input
actspace=`grep "Active_space" ../Infos.dat | awk '{ print $2 }'`
if [[ $actspace == "custom" ]]; then
   actelec=`grep "Active_electrons" ../Infos.dat | awk '{ print $2 }'`
   actorb=`grep "Active_orbitals" ../Infos.dat | awk '{ print $2 }'`
   inactorb=`grep "Inactive_orbitals" ../Infos.dat | awk '{ print $2 }'`

   sed -i "s|nActEl=10 0 0|nActEl=$actelec 0 0|" CASPT2_ipea_025/${Project}_CASPT2_025.input
   sed -i "s|Ras2=10|Ras2=$actorb|" CASPT2_ipea_025/${Project}_CASPT2_025.input
   sed -i "s|Inactive=62|Inactive=$inactorb|" CASPT2_ipea_025/${Project}_CASPT2_025.input

   sed -i "s|nActEl=10 0 0|nActEl=$actelec 0 0|" CASPT2_ipea_0/${Project}_CASPT2_0.input
   sed -i "s|Ras2=10|Ras2=$actorb|" CASPT2_ipea_0/${Project}_CASPT2_0.input
   sed -i "s|Inactive=62|Inactive=$inactorb|" CASPT2_ipea_0/${Project}_CASPT2_0.input
fi

cp $Project_new/$Project_new.key CASPT2_ipea_025/${Project}_CASPT2_025.key
cp $Project_new/$Project_new.key CASPT2_ipea_0/${Project}_CASPT2_0.key

cp $Project_new/$Project_new.xyz CASPT2_ipea_025/${Project}_CASPT2_025.xyz
cp $Project_new/$Project_new.xyz CASPT2_ipea_0/${Project}_CASPT2_0.xyz

cp $Project_new/$prm.prm CASPT2_ipea_025
cp $Project_new/$prm.prm CASPT2_ipea_0

cp $Project_new/$Project_new.JobIph CASPT2_ipea_025/${Project}_CASPT2_025.JobIph
cp $Project_new/$Project_new.JobIph CASPT2_ipea_0/${Project}_CASPT2_0.JobIph

#slurm
cp $Project_new/molcas-job.sh CASPT2_ipea_025
cp $Project_new/molcas-job.sh CASPT2_ipea_0
sed -i "s/NOMEPROGETTO/$Project/" CASPT2_ipea_025/molcas-job.sh
sed -i "s/NOMEPROGETTO/$Project/" CASPT2_ipea_0/molcas-job.sh
sed -i "s/MEMTOT/23000/" CASPT2_ipea_025/molcas-job.sh
sed -i "s/MEMTOT/23000/" CASPT2_ipea_0/molcas-job.sh
sed -i "s/MEMORIA/20000/" CASPT2_ipea_025/molcas-job.sh
sed -i "s/MEMORIA/20000/" CASPT2_ipea_0/molcas-job.sh
sed -i "s/walltime=140/walltime=230/" CASPT2_ipea_025/molcas-job.sh
sed -i "s/walltime=140/walltime=230/" CASPT2_ipea_0/molcas-job.sh

sed -i "s/export Project=.*/export Project=${Project}_CASPT2_025/g" CASPT2_ipea_025/molcas-job.sh
sed -i "s/export Project=.*/export Project=${Project}_CASPT2_0/g" CASPT2_ipea_0/molcas-job.sh

#
# Submiting the job
#
cd CASPT2_ipea_025
#TMPFILE=`mktemp -d /scratch/photon_CASPT2_XXXXXX`
#../../update_infos.sh "CASPT2_iRODS" $TMPFILE ../../Infos.dat
#sed -i "s|TEMPFOLDER|$TMPFILE|g" molcas-job.sh
#cp -r * $TMPFILE
#current=$PWD
#cd $TMPFILE
sbatch molcas-job.sh
#cd $current
cd ..

#cd CASPT2_ipea_0
#sbatch molcas-job.sh

#cd ../../
#cp $templatedir/ASEC/Energies_CAV_TINKER.sh .

