#!/bin/bash
#
# Reading information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
multichain=`grep "MultChain" Infos.dat | awk '{ print $2 }'`
step=`grep "Step" Infos.dat | awk '{ print $2 }'`
charge=`grep "Init_Charge" Infos.dat | awk '{ print $2 }'`
relaxpr=`grep "Relax_protein" Infos.dat | awk '{ print $2 }'`
moldy=`grep "MD_ensemble" Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`

cd Dynamic

#
# Colecting the files needed to run the NPT molecular Dynamics
#
if [[ $step -eq 0 ]]; then
   if [[ $moldy == "NVT" ]]; then
      mkdir Sim_NPT     
      cp $templatedir/ASEC/dynamic_sol_NPT.mdp Sim_NPT
      cp Minimization/final-${Project}_box_sol.gro Sim_NPT/${Project}_box_sol.gro
      cp Minimization/*.itp Sim_NPT
      cp Minimization/${Project}_box_sol.ndx Sim_NPT
      cp Minimization/${Project}_box_sol.top Sim_NPT
      cp -r Minimization/$amber.ff Sim_NPT
      cp Minimization/residuetypes.dat Sim_NPT
      cd Sim_NPT
   else
      cp $templatedir/ASEC/dynamic_sol_NPT.mdp .
      cp Minimization/final-${Project}_box_sol.gro ${Project}_box_sol.gro
      cp Minimization/*.itp .
      cp Minimization/${Project}_box_sol.ndx .
      cp Minimization/${Project}_box_sol.top .
   fi
fi
if [[ $relaxpr == y ]]; then
   sed -i "s/;freezegrps = GroupDyna/freezegrps = GroupDyna/g" dynamic_sol_NPT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NPT.mdp
else
   sed -i "s/;freezegrps = non-Water/freezegrps = non-Water/g" dynamic_sol_NPT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NPT.mdp
fi
#else
#   cd Dynamic
#   sed -i "s/;freezegrps = GroupDyna/freezegrps = GroupDyna/g" dynamic_sol_NPT.mdp
#   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NPT.mdp
#fi

#
# Defining parameters for the MD
#
echo ""  
echo " What is the PRODUCTION TEMPERATURE of the NPT simulation? (Kelvin)"
echo ""
read tempmd
echo ""
echo " Do you want to heat the system before the MD production run? (y/n)"
echo
read risposta
if [[ $risposta == y ]]; then
   echo ""
   echo " How long is the HEATING PHASE? (ps)"
   echo ""
   read timeheat
   echo ""
   echo " How long is the EQUILIBRATION PHASE? (ps)"
   echo ""
   read timequi
   echo ""  
else
   timeheat=0
   timequi=0
fi
echo " How long is the production phase? (ps)"
read timeprod
echo ""  

if [[ $risposta == y ]]; then
   numsteps=$(($timeheat+$timequi+$timeprod))
   sed -i "s/TIME1/$timeheat/" dynamic_sol_NPT.mdp
   sed -i "s/TEMP1/$tempmd/g" dynamic_sol_NPT.mdp
else
   numsteps=$timeprod
   sed -i "s/annealing/;annealing/" dynamic_sol_NPT.mdp
   sed -i "s/;gen_vel/gen_vel/" dynamic_sol_NPT.mdp
   sed -i "s/;gen_temp/gen_temp/" dynamic_sol_NPT.mdp
#   sed -i "s/;gen_temp/gen_temp/" dynamic.mdp
   sed -i "s/ref_t = 0/;ref_t = 0/" dynamic_sol_NPT.mdp
   sed -i "s/;ref_t = TEMP1/ref_t = TEMP1/" dynamic_sol_NPT.mdp
   sed -i "s/TEMP1/$tempmd/g" dynamic_sol_NPT.mdp
fi
numsteps=$(($numsteps*1000))
sed -i "s/PASSI/$numsteps/" dynamic_sol_NPT.mdp

#
# Run in the CPUs or GPUs
#
gpu="b"
while [[ $gpu != "y" && $gpu != "n" ]]; do
   echo ""
   echo ""
   echo " Do you want to use the GPUs to compute the dynamics? (y/n)"
   echo ""
   read gpu
done

if [[ $gpu == y"" ]]; then
   cp $templatedir/gromacs.slurm_GPU.sh gromacs.sh
else
   cp $templatedir/gromacs.slurm.sh gromacs.sh
fi

$gropath/gmx grompp -maxwarn 2 -f dynamic_sol_NPT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr

sed -i "s|NOMEPROGETTO|${Project}_box_sol|" gromacs.sh
sed -i "s|NOMEDIRETTORI|$PWD|" gromacs.sh
sed -i "s|GROPATH|$gropath|" gromacs.sh

sbatch gromacs.sh

#
# Message to user
#

if [[ $moldy == "NVT" ]]; then
   cd ../../
   cp $templatedir/ASEC/MD_NVT.sh .
   cp $templatedir/Analysis_MD.sh .
   echo ""
   echo "**************************************************************"
   echo ""
   echo " Wait for the NPT molecular dynamics to end then run MD_NVT.sh"
   echo ""
   echo "**************************************************************"
   ./update_infos.sh "Heat_NPT" $timeheat Infos.dat
   ./update_infos.sh "Equi_NPT" $timequi Infos.dat
   ./update_infos.sh "Prod_NPT" $timeprod Infos.dat
   ./update_infos.sh "Next_script" "MD_NVT.sh" Infos.dat
else
   cd ../
   cp $templatedir/ASEC/MD_ASEC.sh .
   cp $templatedir/Analysis_MD.sh .
   echo ""
   echo "***************************************************************"
   echo ""
   echo " Wait for the NPT molecular dynamics to end then run MD_ASEC.sh"
   echo ""
   echo "***************************************************************"
   ./update_infos.sh "HeatMD" $timeheat Infos.dat
   ./update_infos.sh "EquiMD" $timequi Infos.dat
   ./update_infos.sh "ProdMD" $timeprod Infos.dat
   ./update_infos.sh "Next_script" "MD_ASEC.sh" Infos.dat
fi

