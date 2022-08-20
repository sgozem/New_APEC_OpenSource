#!/bin/bash
#
# This script generates a PDB file of the final structure
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" ../Infos.dat | awk '{ print $2 }'`

###############################################
###############################################
#    CAVITY ENERGY
###############################################
###############################################

mkdir Total_CASPT2
cp Total_CASSCF/Rerun_1/md.log Total_CASPT2/md_1.log
cp Total_CASSCF/Rerun_2/md.log Total_CASPT2/md_2.log
#cp Total_CASSCF/Rerun_3/md.log Total_CASPT2/md_3.log
#cp ../Dynamic/output/md.log Total_CASPT2

cp CASPT2_ipea_025/${Project}_CASPT2_025.out Total_CASPT2/molcas.out 
cd Total_CASPT2
cp $templatedir/ASEC/Energies_CAV_CASPT2.f .
#cp ../Energies_CAV.f .

answer=0
while  [[ $answer -ne 1 && $answer -ne 2 && $answer -ne 3 ]]; do
   echo ""
   echo ""
   echo " Which CASPT2 root do you want to compute  "
   echo ""
   echo ""
   read answer
done

if [[ $answer -eq 1 ]]; then
   echo " Root number 1 will be computed "
fi
if [[ $answer -eq 2 ]]; then
   sed -i "s/croot_2/\ /g" Energies_CAV_CASPT2.f
   echo " Root number 2 will be computed "
fi
if [[ $answer -eq 3 ]]; then
   sed -i "s/croot_2/\ /g" Energies_CAV_CASPT2.f
   sed -i "s/croot_3/\ /g" Energies_CAV_CASPT2.f
   echo " Root number 3 will be computed "
fi

gfortran Energies_CAV_CASPT2.f -o Energies_CAV_CASPT2.x
./Energies_CAV_CASPT2.x
rm Energies_CAV_CASPT2.x Energies_CAV_CASPT2.f

cd ../

###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

