#!/bin/bash
#OAR -n NOMEPROGETTO
#OAR -l nodes=1,walltime=5:00:00
#--------------------------------------------------------------------#
module load gcc/4.8.3
#--------------------------------------------------------------------#
Project=NOMEPROGETTO
WorkDir=/tmp/$USER/$Project.$OAR_JOB_ID
InpDir=NOMEDIRETTORI
OutDir=$InpDir/output
#--------------------------------------------------------------------#
mkdir -p $WorkDir
mkdir $OutDir
cp $InpDir/* $WorkDir
cd $WorkDir
GROPATH/mdrun -s $Project.tpr -o $Project.trr -x $Project.xtc -c final-$Project.gro
cp * $OutDir
cd $InpDir
rm $WorkDir -r
