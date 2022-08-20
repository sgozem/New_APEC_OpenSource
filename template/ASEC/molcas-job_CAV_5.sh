#!/bin/bash
#PBS -N NOMEPROGETTO
#PBS -S /bin/sh
#PBS -r n 
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2000mb
#PBS -V
#PBS -A PAA0009
#--------------------------------------------------------------#
NPROCS=`wc -l < $PBS_NODEFILE`
cd $PBS_O_WORKDIR
#--------------------------------------------------------------#
# Molcas settings 
#--------------------------------------------------------------#
module load intel-compilers-10.0.023
export MOLCAS="/nfs/03/bgs0361/bin/7.8.dev"
export MOLCASMEM=2000mb
export MOLCAS_MOLDEN=ON
export MOLCAS_PRINT=normal
export TINKER="/nfs/03/bgs0361/bin/7.8.dev/tinker/bin_qmmm"
#--------------------------------------------------------------#
#  Change the Project!!!
#--------------------------------------------------------------#
echo $HOSTNAME > $InpDir/nodename
echo $JOBID > $InpDir/jobid
export WorkDir=/tmp/$Project.$PBS_JOBID
mkdir -p $WorkDir
export InpDir=$PBS_O_WORKDIR
cd $WorkDir
for i in {401..500}
do
  export Project=Tinker_cav_$i
#--------------------------------------------------------------#
# Copy of the files - obsolete
#--------------------------------------------------------------#
#cp $InpDir/$Project.xyz $WorkDir/$Project.xyz
#cp $InpDir/$Project.key $WorkDir/$Project.key
#cp $InpDir/*.prm $WorkDir/
#--------------------------------------------------------------#
# Start job
#--------------------------------------------------------------#
  /nfs/03/bgs0361/bin/dowser/bin/molcas $InpDir/$Project.input >$InpDir/$Project.out 2>$InpDir/$Project.err
  rm $WorkDir/*
done

