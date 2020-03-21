#!/bin/bash
#PBS -N NOMEPROGETTO
#PBS -S /bin/sh
#PBS -r n 
#PBS -l walltime=hh:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=MEMTOTMB
#PBS -V
#PBS -A PAA0009
#PBS -m ae
#--------------------------------------------------------------#
NPROCS=`wc -l < $PBS_NODEFILE`
cd $PBS_O_WORKDIR
#--------------------------------------------------------------#
# Molcas settings 
#--------------------------------------------------------------#
module load intel-compilers-10.0.023
export MOLCAS="/users/PCS0202/bgs0361/bin/7.8.dev"
export MOLCASMEM=MEMORIAMB
export MOLCAS_MOLDEN=ON
export MOLCAS_PRINT=normal
export TINKER="/users/PCS0202/bgs0361/bin/7.8.dev/tinker/bin_qmmm"
#--------------------------------------------------------------#
#  Change the Project!!!
#--------------------------------------------------------------#
export Project=$PBS_JOBNAME
export WorkDir=/tmp/$Project.$PBS_JOBID
mkdir -p $WorkDir
export InpDir=$PBS_O_WORKDIR
echo $HOSTNAME > $InpDir/nodename
echo $JOBID > $InpDir/jobid
#--------------------------------------------------------------#
# Copy of the files - obsolete
#--------------------------------------------------------------#
#cp $InpDir/$Project.xyz $WorkDir/$Project.xyz
#cp $InpDir/$Project.key $WorkDir/$Project.key
#cp $InpDir/*.prm $WorkDir/
#--------------------------------------------------------------#
# Start job
#--------------------------------------------------------------#
cd $WorkDir
/users/PCS0202/bgs0361/bin/dowser/bin/molcas $InpDir/$Project.input >$InpDir/$Project.out 2>$InpDir/$Project.err
