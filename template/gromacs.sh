#!/bin/bash
#PBS -N NOMEPROGETTO
#PBS -S /bin/sh
#PBS -r n 
#PBS -l walltime=60:00:00
#PBS -l mem=40GB
#PBS -l nodes=1:ppn=16
#PBS -A PAA0009
#PBS -m ae
#--------------------------------------------------------------#
NPROCS=`wc -l < $PBS_NODEFILE`
module load gromacs/4.6.3
module load modules/au2014
cd $PBS_O_WORKDIR
export Project=$PBS_JOBNAME
export WorkDir=$TMPDIR/$Project.$PBS_JOBID
export InpDir=NOMEDIRETTORI
export outdir=NOMEDIRETTORI/output
echo $HOSTNAME > $InpDir/nodename
echo $JOBID > $InpDir/jobid
mkdir $outdir
mkdir -p $WorkDir
#--------------------------------------------------------------#
# Start job
#--------------------------------------------------------------#
cp $InpDir/* $WorkDir
cd $WorkDir
mpiexec GROPATH/mdrun_mpi -s $Project.tpr -o $Project.trr -x $Project.xtc -c final-$Project.gro
cp * $outdir/


