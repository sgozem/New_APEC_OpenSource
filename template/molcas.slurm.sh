#!/bin/bash
#SBATCH -J NOMEPROGETTO
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t hh:00:00
#SBATCH --mem=MEMTOTMB 
#SBATCH -p qCPU120
#SBATCH -A CHEM9C4
##SBATCH --exclude photongpu01,photongpu02
#--------------------------------------------------------------#
#cd $SLURM_SUBMIT_DIR
#--------------------------------------------------------------#
# Molcas settings 
#--------------------------------------------------------------#
export MOLCAS="/userapp/APEC_KABIR/build"
export MOLCAS_MEM=MEMORIAMB
export MOLCAS_MOLDEN=ON
export MOLCAS_PRINT=normal
export TINKER="/userapp/APEC_KABIR/build/tinker/bin"
export WorkDir=/scratch/$SLURM_JOB_ID
export InpDir=$PWD
#export outdir=$PWD/output
export Project=$SLURM_JOB_NAME
#export WorkDir=/runjobs/RS10237/$SLURM_JOB_ID
#--------------------------------------------------------------#
#  Change the Project!!!
#--------------------------------------------------------------#
mkdir -p $WorkDir
#mkdir $outdir
#export InpDir=$SLURM_SUBMIT_DIR
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
cp $InpDir/* $WorkDir
cd $WorkDir
/userapp/APEC_KABIR/bin/pymolcas $WorkDir/$Project.input >$WorkDir/$Project.out 2>$WorkDir/$Project.err
#cp $WorkDir/* $outdir
cp $WorkDir/* $InpDir
rm $InpDir/*.OrdInt
rm -r $WorkDir

