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

config=1000

cd Total_TINKER
if [ -f CAV_energies ]; then
   rm CAV_energies
fi
cd Tinker_calc
for i in $(eval echo "{1..$config}")
do
   grep " eb   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ea   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eba  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eub  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eaa  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eopb = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eopd = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eid  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eit  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " et   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ebt  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ett  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ev   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ec   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ecd  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ed   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " em   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ep   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " er   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " es   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " elf  = " Tinker_cav_$i.out >> ../CAV_energies
   grep " eg   = " Tinker_cav_$i.out >> ../CAV_energies
   grep " ex   = " Tinker_cav_$i.out >> ../CAV_energies
done
cd ..
cp ../CASPT2_ipea_025/${Project}_CASPT2_025.out molcas.out 
cp $templatedir/ASEC/Total_TINKER.f .

iterms=`grep -w -c " eb \| ea \| eba \| eub \| eaa \| eopb \| eopd \| eid \| eit \| et \| ebt \| ett \| ev \| ec \| ecd \| ed \| em \| ep \| er \| es \| elf \| eg \| ex " Tinker_calc/Tinker_cav_1.out`
iterms2=`grep -w -c " eb \| ea \| eba \| eub \| eaa \| eopb \| eopd \| eid \| eit \| et \| ebt \| ett \| ev \| ec \| ecd \| ed \| em \| ep \| er \| es \| elf \| eg \| ex " molcas.out`

sed -i "s/terminos2/$iterms2/g" Total_TINKER.f
sed -i "s/terminos/$iterms/g" Total_TINKER.f

sed -i "s/configu/$config/g" Total_TINKER.f

lines=`wc -l CAV_energies | awk '{ print $1 }'`
sed -i "s/lineas/$lines/g" Total_TINKER.f

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
   sed -i "s/croot_2/\ /g" Total_TINKER.f
   echo " Root number 2 will be computed "
fi
if [[ $answer -eq 3 ]]; then
   sed -i "s/croot_2/\ /g" Total_TINKER.f
   sed -i "s/croot_3/\ /g" Total_TINKER.f
   echo " Root number 3 will be computed "
fi

gfortran Total_TINKER.f -o Total_TINKER.x
./Total_TINKER.x
#rm Total_TINKER.x Total_TINKER.f

cd ../

###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

