#!/bin/bash
#OAR -n NOMEPROGETTO
#OAR -l nodes=1,walltime=hh:00:00
#--------------------------------------------------------------------#
module load gcc/4.8.3
#--------------------------------------------------------------------#
export MOLCAS="/nfs/03/bgs0361/bin/7.8.dev"
export MOLCAS_MEM=MEMORIAMB
export MOLCAS_MOLDEN=ON
export MOLCAS_PRINT=normal
export TINKER="/nfs/03/bgs0361/bin/7.8.dev/tinker/bin_qmmm"
export Project=NOMEPROGETTO
export WorkDir=/tmp/$USER/$Project.$OAR_JOB_ID
export InpDir=NOMEDIRETTORI
#--------------------------------------------------------------------#
mkdir -p $WorkDir
cd $WorkDir
/nfs/03/bgs0361/bin/dowser/bin/molcas $InpDir/$Project.input >$InpDir/$Project.out 2>$InpDir/$Project.err
cd $InpDir
rm $WorkDir -r
