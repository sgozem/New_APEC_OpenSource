#!/bin/bash
#SBATCH -J NOMEPROGETTO
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t hh:00:00
#SBATCH --mem=MEMTOTMB 
#SBATCH -p qPHO
#SBATCH --exclude photongpu01,photongpu02
#--------------------------------------------------------------#
cd $SLURM_SUBMIT_DIR
#--------------------------------------------------------------#
# Molcas settings 
#--------------------------------------------------------------#
export MOLCAS="/home/users/yorozcogonzalez/bin/build"
export MOLCAS_MEM=MEMORIAMB
export MOLCAS_MOLDEN=ON
export MOLCAS_PRINT=normal
export TINKER="/home/users/yorozcogonzalez/bin/build/tinker/bin"
#--------------------------------------------------------------#
#  Change the Project!!!
#--------------------------------------------------------------#
export Project=$SLURM_JOB_NAME
export WorkDir=/runjobs/RS10237/$SLURM_JOB_ID
mkdir -p $WorkDir
export InpDir=$SLURM_SUBMIT_DIR
echo $SLURM_JOB_NODELIST > $InpDir/nodename
echo $SLURM_JOB_ID > $InpDir/jobid
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
/home/users/yorozcogonzalez/bin/pymolcas $InpDir/$Project.input >$InpDir/$Project.out 2>$InpDir/$Project.err
cp $WorkDir/$Project.RasOrb $InpDir
cp $WorkDir/$Project.*.molden $InpDir
cp $WorkDir/$Project.*.Tinker.log $InpDir
rm -r $WorkDir

