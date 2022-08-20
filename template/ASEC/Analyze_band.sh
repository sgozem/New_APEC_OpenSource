#!/bin/bash
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`

cd Band_width

confs=100

echo "  Root 1          Root 2          Root 3         Absorption" > Energies
for i in $(eval echo "{1..$confs}"); do
   S0=`grep "::    CASPT2 Root  1     Total energy:" ${Project}_conf_${i}.out | awk '{ print $7 }'`
   S1=`grep "::    CASPT2 Root  2     Total energy:" ${Project}_conf_${i}.out | awk '{ print $7 }'`
   S2=`grep "::    CASPT2 Root  3     Total energy:" ${Project}_conf_${i}.out | awk '{ print $7 }'`

   abs=$(echo "($S1 - $S0)*627.5091809" | bc)
   echo " $S0   $S1   $S2   $abs" >> Energies
done

cp Energies ../

