#!/bin/bash
#
# This script generates a PDB file of the final structure
# Retrieving information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
seeds=`grep "Parallel_MD" Infos.dat | awk '{ print $2 }'`
heat=`grep "HeatMD" Infos.dat | awk '{ print $2 }'`
equi=`grep "EquiMD" Infos.dat | awk '{ print $2 }'`
prod=`grep "ProdMD" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`
relaxpr=`grep "Relax_protein" Infos.dat | awk '{ print $2 }'`

###############################################
###############################################
#    CAVITY ENERGY
###############################################
###############################################

echo ""
echo "***********************************************************"
echo ""
echo " The total MM energy of the system will be computed by"
echo " re-runing the MD with ZERO paremeters of the flavin"
echo ""
echo "***********************************************************"
echo ""

mkdir calculations/Total_QMMM_Energy
cd calculations/Total_QMMM_Energy

for k in $(eval echo "{1..$seeds}")
do
   echo ""
   echo " Wait for seed $k ..."
   mkdir Rerun_$k
   cd Rerun_$k
 
   cp ../../../Dynamic/seed_$k/*.itp .
   cp ../../../Dynamic/seed_$k/*.top .
   cp ../../../Dynamic/seed_$k/dynamic_sol_NVT.mdp .
   cp ../../../Dynamic/seed_$k/${Project}_box_sol.gro .
   cp ../../../Dynamic/seed_$k/${Project}_box_sol.ndx .
   cp ../../../Dynamic/seed_$k/output/${Project}_box_sol.trr .
   cp -r ../../../Dynamic/seed_$k/amber99sb.ff .

   sed -i "s/annealing/;annealing/" dynamic_sol_NVT.mdp
   sed -i "s/annealing_npoints/;annealing_npoints/" dynamic_sol_NVT.mdp
   sed -i "s/annealing_time/;annealing_time/" dynamic_sol_NVT.mdp
   sed -i "s/annealing_temp/;annealing_temp/" dynamic_sol_NVT.mdp
   sed -i "s/ref_t = 0/;ref_t = 0/" dynamic_sol_NVT.mdp
   sed -i "s/;ref_t = 300/ref_t = 300/" dynamic_sol_NVT.mdp

   sed -i 's/ C\* / XX /g' ${Project}_Other_chain_A2.itp
   sed -i 's/ N\* / YY /g' ${Project}_Other_chain_A2.itp

   grep "     1    CHR   " ${Project}_Other_chain_A2.itp > temp1

   sed -i 's/nstlog                  = 100/nstlog                  = 500/g' dynamic_sol_NVT.mdp

   lines=`wc -l temp1 | awk '{ print $1 }'`
   for i in $(eval echo "{1..$lines}")
   do
      head -n $i temp1 | tail -n1 > temp2
      line1=`head -n1 temp2 | awk '{ print $0 }'`
      carga=`head -n1 temp2 | awk '{ print $7 }'`
      atomt=`head -n1 temp2 | awk '{ print $2 }'`
      sed -i "s/$carga/ 0.0000/" temp2
      sed -i "s/  $atomt/CH${atomt}/" temp2
      line2=`head -n1 temp2 | awk '{ print $0 }'`
      sed -i "s/$line1/$line2/" ${Project}_Other_chain_A2.itp
      sed -i "s/ CHXX / CHC\* /g" ${Project}_Other_chain_A2.itp
      sed -i "s/ CHYY / CHN\* /g" ${Project}_Other_chain_A2.itp
   done

   cp $templatedir/amber99sb.ff/CHffbonded.itp amber99sb.ff/ffbonded.itp
   cp $templatedir/amber99sb.ff/CHffnonbonded.itp amber99sb.ff/ffnonbonded.itp

   cp $templatedir/gromacs.slurm.rerun.sh gromacs.sh
   sed -i "s|NOMEPROGETTO|${Project}_box_sol|" gromacs.sh
   sed -i "s|NOMEDIRETTORI|$PWD|" gromacs.sh
   sed -i "s|GROPATH|$gropath|" gromacs.sh

   sed -i "s/SBATCH -t 23:59:00/SBATCH -t 10:00:00/" gromacs.sh

   $gropath/gmx grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr
   sbatch gromacs.sh

#   $gropath/grompp -maxwarn 2 -f dynamic_sol_NVT.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr
#   $gropath/mdrun -s ${Project}_box_sol.tpr -rerun ${Project}_box_sol.trr
#   cp md.log ../md_rerun_$k.log

   cd ..
cp $templatedir/ASEC/Total_QMMM_Energy_results.sh ../
done

echo ""
echo "**********************************************************"
echo ""
echo " The re-run calculations were submited. When this is done"
echo " run Total_QMMM_next.sh"
echo ""
echo "**********************************************************"
echo ""

