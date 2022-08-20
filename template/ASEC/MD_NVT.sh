#!/bin/bash
#
# Reading project name, parm file and templatedir from Infos.dat
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
tmpdir=`grep "tempdir" Infos.dat | awk '{ print $2 }'`

cd Dynamic
if [[ $step -eq 0 ]]; then
   if [[ -d Sim_NPT/output ]]; then
      echo ""
      echo " *********************************************************************"
      echo "                      Warning!"
      echo ""
      echo " MD Sim_NPT/output directory already excist. We are goint to use it..."
      echo ""
      echo " *********************************************************************"
      echo ""
   else
      echo ""
      echo " Collecting data from the nodes, please wait ..."
      echo ""
      dir=`basename $tmpdir`
      iget -r /arctic/projects/CHEM9C4/$USER/$dir Sim_NPT
      if [[ -f Sim_NPT/$dir/md.log ]]; then
        mv Sim_NPT/$dir Sim_NPT/output 
        irm -r /arctic/projects/CHEM9C4/$USER/$dir
      else
         echo ""
         echo "************************************************************************"
         echo ""
         echo " It seems that the NPT MD is still running or it did not finish properly"
         echo ""
         echo "************************************************************************"
         echo ""
         exit 0
      fi
   fi
   cp $templatedir/ASEC/dynamic_sol_NVT.mdp .
   cp Sim_NPT/output/final-${Project}_box_sol.gro ${Project}_box_sol.gro
   cp Sim_NPT/*.itp .
   cp Sim_NPT/${Project}_box_sol.ndx .
   cp Sim_NPT/${Project}_box_sol.top .
   cp Sim_NPT/residuetypes.dat .
   cp -r Sim_NPT/$amber.ff .
fi

if [[ $relaxpr == y ]]; then
   sed -i "s/;freezegrps = GroupDyna/freezegrps = GroupDyna/g" dynamic_sol_NVT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NVT.mdp
else
   sed -i "s/;freezegrps = non-Water/freezegrps = non-Water/g" dynamic_sol_NVT.mdp
   sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" dynamic_sol_NVT.mdp
fi
echo ""  
echo " What is the PRODUCTION TEMPERATURE of the NVT simulation? (Kelvin). Normally 300."
echo ""
read tempmd

if [[ $step -eq 0 ]]; then
   heating="n"
else
   #echo ""
   #echo " Do you want to heat the system before the MD production run? (y/n)"
   #echo
   #read risposta
   heating="y"
fi
if [[ $heating == y ]]; then
   echo ""
   echo " The system will to be heated from 0k to the defined temperature"
   echo " because at this point we do not have velocities from previous calculations."
   echo ""
   echo " How long is the HEATING PHASE? (ps). Normally 300."
   echo ""
   read timeheat
   echo ""
   echo " How long is the EQUILIBRATION PHASE? (ps). Normally 4700."
   echo ""
   read timequi
   echo ""
else
   echo ""
   echo " Here we do not need to heat the system because we have"
   echo " the velocities from the NPT Molecular dynamics."
   echo ""
   echo " How long is the EQUILIBRATION PHASE? (ps). Normally 5000"
   echo ""
   read timequi
   timeheat=0
   #timequi=0
fi
echo ""
echo " How long is the PRODUCTION PHASE? (ps). Normally 5000 in 2 seeds if in Step 0"
echo ""
read timeprod
echo ""

parallelize=r
while [[ $parallelize != y && $parallelize != n ]]; do
   echo ""
   echo " Do you want to parallelize the production phase of the MD? (y/n)"
   echo ""
   echo " Take into account how many GPU nodes are available if you want to"
   echo " use the GPUs to run the Molecular Dynamics."
   echo ""
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

if [[ $heating == y ]]; then
   numsteps=$(($timeheat+$timequi+$timeprod))
   sed -i "s/TIME1/$timeheat/" dynamic_sol_NVT.mdp
   sed -i "s/TEMP1/$tempmd/g" dynamic_sol_NVT.mdp
else
   numsteps=$(($timequi+$timeprod))
   sed -i "s/annealing/;annealing/" dynamic_sol_NVT.mdp
#   sed -i "s/;gen_vel/gen_vel/" dynamic_sol_NVT.mdp    # it is goint to read vel from last .gro
#   sed -i "s/;gen_temp/gen_temp/" dynamic_sol_NVT.mdp
#   sed -i "s/;gen_temp/gen_temp/" dynamic.mdp
   sed -i "s/ref_t = 0/;ref_t = 0/" dynamic_sol_NVT.mdp
   sed -i "s/;ref_t = TEMP1/ref_t = TEMP1/" dynamic_sol_NVT.mdp
   sed -i "s/TEMP1/$tempmd/g" dynamic_sol_NVT.mdp
fi
numsteps=$(($numsteps*1000))
sed -i "s/PASSI/$numsteps/" dynamic_sol_NVT.mdp

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

      TMPFILE=`mktemp -d /scratch/photon_seed_${i}_XXXXXX`
      ../../update_infos.sh "MD_seed_$i" $TMPFILE ../../Infos.dat
      sed -i "s|TEMPFOLDER|$TMPFILE|" gromacs.sh
      cp -r * $TMPFILE
      current=$PWD
      cd $TMPFILE
      sbatch gromacs.sh
      cd $current

      cd ..
      done
   else # no parallel
      $gropath/gmx grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr

      sed -i "s|NOMEPROGETTO|${Project}_box_sol|" gromacs.sh
      sed -i "s|NOMEDIRETTORI|$PWD|" gromacs.sh
      sed -i "s|GROPATH|$gropath|" gromacs.sh

      sed -i "s/SBATCH -t 23:59:00/SBATCH -t 47:59:00/" gromacs.sh

      TMPFILE=`mktemp -d /scratch/photon_XXXXXXXX`
      ../update_infos.sh "tempdir" $TMPFILE ../Infos.dat
      sed -i "s|TEMPFOLDER|$TMPFILE|" gromacs.sh
      cp -r * $TMPFILE
      current=$PWD
      cd $TMPFILE
      sbatch gromacs.sh
      cd $current

fi

cd ../
cp $templatedir/ASEC/MD_ASEC.sh .
./update_infos.sh "Next_script" "MD_ASEC.sh" Infos.dat

cp $templatedir/Analysis_MD.sh .
echo ""
echo "***************************************************************"
echo ""
echo " Wait for the NVT molecular dynamics to end then run MD_ASEC.sh"
echo ""
echo "***************************************************************"
