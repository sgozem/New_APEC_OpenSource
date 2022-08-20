#!/bin/bash
#SBATCH -J NOMEPROGETTO
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 23:59:00
#SBATCH --mem=45G 
#SBATCH -p qPHOGPU
#SBATCH --gres=gpu:p100:1
#--------------------------------------------------------------#
#NPROCS=`wc -l < $SLURM_JOB_NODELIST`
cd $SLURM_SUBMIT_DIR

module load ComputationalChemistry/Gromacs2019

export Project=$SLURM_JOB_NAME
export WorkDir=/runjobs/RS10237/$SLURM_JOB_ID
export InpDir=NOMEDIRETTORI
export outdir=NOMEDIRETTORI
echo $SLURM_JOB_NODELIST > $InpDir/nodename
echo $SLURM_JOB_ID > $InpDir/jobid
#mkdir $outdir
mkdir -p $WorkDir
#--------------------------------------------------------------#
# Start job
#--------------------------------------------------------------#
cp -r $InpDir/* $WorkDir
cd $WorkDir
#/home/users/yorozcogonzalez/bin/gromacs/bin
#/apps/Computation_Chemistry/gromacs-2019/Debug_Build-Skylake-GPU/bin/gmx grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr
/apps/Computation_Chemistry/gromacs-2019/Debug_Build-Skylake-GPU/bin/gmx mdrun -nt 16 -s $Project.tpr -rerun $Project.trr
cp $WorkDir/* $outdir/

rm -r $WorkDir

