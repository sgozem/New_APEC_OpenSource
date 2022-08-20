#!/bin/bash
#
# This script generates a PDB file of the final structure
# Retrieving information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`

  cd MD_analysis
  mkdir AutoCorrelation
  cd AutoCorrelation
  mkdir new_Dynamic
  cd new_Dynamic

cp -r ../../../Dynamic/amber94.ff .
cp ../../../Dynamic/*.itp .
cp ../../../Dynamic/$Project.gro .
cp ../../../Dynamic/$Project.top .
cp ../../../Dynamic/$Project.ndx .
cp ../../../Dynamic/posre* .
cp ../../../Dynamic/residuetypes.dat .
cp ../../../Dynamic/dynamic.mdp .
cp ../../../Dynamic/gromacs.sh .

sed -i "s/nsteps = .*/nsteps = 400000/g" dynamic.mdp
sed -i "s/nstxout = .*/nstxout = 10/g" dynamic.mdp
sed -i "s/nstvout = .*/nstvout = 10/g" dynamic.mdp
sed -i "s/nstlog = .*/nstlog = 10/g" dynamic.mdp
sed -i "s|export inpdir=.*|export inpdir="$(pwd)"|g" gromacs.sh
sed -i "s|export outdir=.*|export outdir="$(pwd)"/output|g" gromacs.sh

$gropath/gmx grompp -f dynamic.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr > mdlog

cp $templatedir/ASEC/Autocorrelation_2.sh ../../../

qsub gromacs.sh

echo ""
echo " Wait for new molecular dynamics to end, then run Autocorrelation_2.sh "
echo ""


###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

