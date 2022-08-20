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


  echo ""
  echo ""
  echo " Processing "
  echo ""
  echo ""

  cd MD_analysis
#  mkdir AutoCorrelation
  cd AutoCorrelation
#  mkdir Rerun


echo ""
echo " Wait for rerun to end, then run Autocorrelation_3.sh "
echo ""

cp Rerun/md.log md_zero.log
cp new_Dynamic/output/md.log .
#cd ../

cp $templatedir/ASEC/Autocorrelation.f .

gfortran Autocorrelation.f -o Autocorrelation.x
./Autocorrelation.x
rm Autocorrelation.x

###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

