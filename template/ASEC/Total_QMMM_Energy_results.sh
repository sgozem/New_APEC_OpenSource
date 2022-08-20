#!/bin/bash
#
# This script generates a PDB file of the final structure
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`
seeds=`grep "Parallel_MD" ../Infos.dat | awk '{ print $2 }'`
heat=`grep "HeatMD" ../Infos.dat | awk '{ print $2 }'`
equi=`grep "EquiMD" ../Infos.dat | awk '{ print $2 }'`
prod=`grep "ProdMD" ../Infos.dat | awk '{ print $2 }'`

###############################################
###############################################
#    CAVITY ENERGY
###############################################
###############################################


echo ""
echo ""
echo " The total QMMM energy of the system will be computed"
echo " The folder Total_QMMM_Energy will be created containing"
echo " results in the file Total_QMMM_Energy.dat"
echo ""
echo ""

cd Total_QMMM_Energy
for k in $(eval echo "{1..$seeds}")
do
   cd Rerun_$k
   cp md.log ../md_rerun_$k.log
   cd ..
done

cp ../CASPT2_ipea_025/${Project}_CASPT2_025.out CASPT2.out
cp $templatedir/ASEC/Total_QMMM_Energy.f .

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
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   echo " Root number 2 will be computed "
fi
if [[ $answer -eq 3 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   echo " Root number 3 will be computed "
fi
if [[ $answer -eq 4 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_4/\ /g" Total_QMMM_Energy.f
   echo " Root number 4 will be computed "
fi
if [[ $answer -eq 5 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_4/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_5/\ /g" Total_QMMM_Energy.f
   echo " Root number 5 will be computed "
fi
if [[ $answer -eq 6 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_4/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_5/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_6/\ /g" Total_QMMM_Energy.f
   echo " Root number 6 will be computed "
fi
if [[ $answer -eq 7 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_4/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_5/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_6/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_7/\ /g" Total_QMMM_Energy.f
   echo " Root number 7 will be computed "
fi
if [[ $answer -eq 8 ]]; then
   sed -i "s/croot_2/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_3/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_4/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_5/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_6/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_7/\ /g" Total_QMMM_Energy.f
   sed -i "s/croot_8/\ /g" Total_QMMM_Energy.f
   echo " Root number 8 will be computed "
fi
#sed -i "s/semilla/$seeds/" Total_QMMM_Energy.f
#total=`grep -c "   Energies (kJ/mol)" md_rerun_1.log | awk '{ print $1 }'`
# total-1 because the initoal structure is also printed
#total=$(($total-1))
#thermal=$((($heat+$equi)*($total)/(prod/seeds+$heat+$equi)))
#sed -i "s/todas/$total/" Total_QMMM_Energy.f
#sed -i "s/termalizacion/$thermal/" Total_QMMM_Energy.f
gfortran Total_QMMM_Energy.f -o Total_QMMM_Energy.x
./Total_QMMM_Energy.x


