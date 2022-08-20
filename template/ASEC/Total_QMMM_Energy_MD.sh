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
Step=`grep "Step" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" ../Infos.dat | awk '{ print $2 }'`
relaxpr=`grep "Relax_protein" ../Infos.dat | awk '{ print $2 }'`

if [ -d ../../Step_${Step}_Total_QMMM ]; then
      echo " Folder \"Step_${Step}_Total_QMMM\" found! Something is wrong ..."
      echo " Terminating ..."
      exit 0
      echo ""
fi

new_folder=Step_${Step}_Total_QMMM
mkdir ../../${new_folder}

cd ..

mkdir ../${new_folder}/Dynamic
cp Infos.dat ../${new_folder}
./update_infos.sh "Step" "Total_QMMM" ../${new_folder}/Infos.dat
cp update_infos.sh $prm.prm template_* ../${new_folder}
cp -r Chromophore ../${new_folder}
cp $templatedir/ASEC/dynamic_sol_NVT.mdp ../${new_folder}/Dynamic
cp $templatedir/gromacs.sh ../${new_folder}/Dynamic
cp -r Dynamic/$amber.ff ../${new_folder}/Dynamic
cp Dynamic/${Project}_box_sol.gro ../${new_folder}/Dynamic
cp Dynamic/${Project}_box_sol.top ../${new_folder}/Dynamic
cp Dynamic/${Project}_box_sol.ndx ../${new_folder}/Dynamic
cp Dynamic/residuetypes.dat ../${new_folder}/Dynamic
cp Dynamic/*.itp ../${new_folder}/Dynamic
cp -r calculations ../${new_folder}

cd ../${new_folder}
cd Dynamic

if [[ $relaxpr == y ]]; then
   sed -i "s/;freezegrps = GroupDyna/freezegrps = GroupDyna/g" dynamic_sol_NVT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NVT.mdp
else
   sed -i "s/;freezegrps = non-Water/freezegrps = non-Water/g" dynamic_sol_NVT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NVT.mdp
fi

tempmd=300
timeheat=300
timequi=4700
timeprod=250000

parallelize=r
while [[ $parallelize != y && $parallelize != n ]]; do
   echo ""
   echo "********************************************************************"
   echo ""
   echo " A 250 000 ps MD will be performed, saving configurations every 5 ps"
   echo ""
   echo " Do you want to parallelize the production phase of the MD? (y/n)"
   echo " (Take into account how many GPU nodes are available if you want to"
   echo " use the GPUs to run the Molecular Dynamics)"
   echo ""
   echo "********************************************************************"
   read parallelize
   echo ""
done

numparallel=1
if [[ $parallelize == y ]]; then
   echo ""
   echo " How many MDs in parallel?"
   echo ""
   read numparallel
   echo ""
fi

../update_infos.sh "HeatMD" $timeheat ../Infos.dat
../update_infos.sh "EquiMD" $timequi ../Infos.dat
../update_infos.sh "ProdMD" $timeprod ../Infos.dat
../update_infos.sh "Parallel_MD" $numparallel ../Infos.dat

if [[ $parallelize == y ]]; then
   timeprod=$(($timeprod/$numparallel))
fi

numsteps=$(($timeheat+$timequi+$timeprod))
sed -i "s/TIME1/$timeheat/" dynamic_sol_NVT.mdp
sed -i "s/TEMP1/$tempmd/g" dynamic_sol_NVT.mdp
numsteps=$(($numsteps*1000))
sed -i "s/PASSI/$numsteps/" dynamic_sol_NVT.mdp
sed -i "s/nstxout                 = 50000/nstxout                 = 5000/" dynamic_sol_NVT.mdp

gpu="b"
while [[ $gpu != "y" && $gpu != "n" ]]; do
   echo ""
   echo ""
   echo " Do you want to use the GPUs to compute the dynamics? (y/n)"
   echo ""
   echo " Take into account how many GPU nodes are available if you"
   echo " to parallelize the production next."
   echo ""
   echo ""
   read gpu
done

if [[ $gpu == y"" ]]; then
   cp $templatedir/gromacs.slurm_GPU.sh gromacs.sh
else
   cp $templatedir/gromacs.slurm.sh gromacs.sh
fi

if [[ $parallelize == y ]]; then
   for i in $(eval echo "{1..$numparallel}")
   do
      mkdir seed_$i
      cp -r $amber.ff seed_$i
      cp ${Project}_box_sol.* *.mdp *.itp *.dat *.sh seed_$i

      cd seed_$i
      sed -i "s/;gen_temp/gen_temp                = $tempmd/" dynamic_sol_NVT.mdp
      sed -i "s/;gen_vel/gen_vel/" dynamic_sol_NVT.mdp
      sed -i "/gen_vel/agen_seed                = $((23456+131*$i))" dynamic_sol_NVT.mdp

      $gropath/gmx grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr

      sed -i "s|NOMEPROGETTO|${Project}_box_sol|" gromacs.sh
      sed -i "s|NOMEDIRETTORI|$PWD|" gromacs.sh
      sed -i "s|GROPATH|$gropath|" gromacs.sh

      sed -i "s/SBATCH -t 23:59:00/SBATCH -t 47:59:00/" gromacs.sh

      sbatch gromacs.sh

      cd ..
      done
   else
      $gropath/gmx grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr

      sed -i "s|NOMEPROGETTO|${Project}_box_sol|" gromacs.sh
      sed -i "s|NOMEDIRETTORI|$PWD|" gromacs.sh
      sed -i "s|GROPATH|$gropath|" gromacs.sh

      sed -i "s/SBATCH -t 23:59:00/SBATCH -t 100:00:00/" gromacs.sh

      sbatch gromacs.sh
fi

cp $templatedir/ASEC/Total_QMMM_Energy_rerun.sh ../

echo "***********************************************************************"
echo ""
echo " A new folder nanmed $new_folder was created and the MD submited"
echo " When this is done do to $new_folder and run Total_QMMM_Energy_rerun.sh"
echo ""
echo "***********************************************************************"

