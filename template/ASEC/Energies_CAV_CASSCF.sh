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
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`

###############################################
###############################################
#    CAVITY ENERGY
###############################################
###############################################


cavi=b
while [[ $cavi != y && $cavi != n ]]; do
      echo " Do you already have the files for calculating the energy of the cavity? (y/n)"
      read cavi
done

if [[ $cavi ==  "n" ]]; then

  echo ""
  echo ""
  echo " The files for calculating the energy of the cavity will be prepared "
  echo " This might take a few minutes"
  echo ""
  echo ""

  mkdir Total_CASSCF
  cd Total_CASSCF
  mkdir Rerun_1
  mkdir Rerun_2
#  mkdir Rerun_3
  cp ../../Dynamic/${Project}_Protein_chain_?.itp .
  cp ../../MD_analysis/list_gro.dat .
  
#######################
# In the following loop will be generated two .itp files for each chain. Number 1
# for calculating the total MM energy of the protein excluding the RET+LYS and 
# number 2 for calculating the energy of the protein excluding the RET+LYS and also 
# all the side chains of the cavity. The diference between this two energies is the 
# energy of the cavity.
#######################

  for k in ${Project}_Protein_chain_?.itp
  do

### ${k%.itp} get the basename of the files .itp, for instance:
### chain_A.itp become chain_A_1.itp and chain_A_2.itp
 
    name1=${k%.itp}_1.itp
    name2=${k%.itp}_2.itp
#    name3=${k%.itp}_3.itp

    cp $k $name1
    cp $k $name2
#    cp $k $name3

    sed -i 's/ C\* / XX /g' $name1
    sed -i 's/ N\* / YY /g' $name1
    sed -i 's/ C\* / XX /g' $name2
    sed -i 's/ N\* / YY /g' $name2
#    sed -i 's/ C\* / XX /g' $name3
#    sed -i 's/ N\* / YY /g' $name3


    grep -w " RET " $name1 | grep -w " CB " > temp
    grep -w " RET " $name1 | grep -w " HB. " >> temp
    grep -w " RET " $name1 | grep -w " CG " >> temp
    grep -w " RET " $name1 | grep -w " HG. " >> temp
    grep -w " RET " $name1 | grep -w " CD " >> temp
    grep -w " RET " $name1 | grep -w " HD. " >> temp
    grep -w " RET " $name1 | grep -w " CE " >> temp
    grep -w " RET " $name1 | grep -w " HE. " >> temp
    grep -w " RET " $name1 | grep -w " NZ " >> temp
    grep -w " RET " $name1 | grep -w " HZ1 " >> temp
    grep -w " C3R \| C2R \| HR " $name1 >> temp

# temp1 save only the RET+LYS atoms
    cp temp temp1

# cont number of lines
    lines=`wc -l list_gro.dat | awk '{ print $1 }'`

#    if [ -f temp ]; then
#       rm temp
#    fi

    for i in $(eval echo "{1..$lines}")
    #for i in {1..5} for testing
    do
         amino=`head -n $i list_gro.dat | tail -n1 | awk '{ print $2 }'`
          atom=`head -n $i list_gro.dat | tail -n1 | awk '{ print $3 }'`
      chainnum=`head -n $i list_gro.dat | tail -n1 | awk '{ print $4 }'`
  
    grep -w " $chainnum    $amino " $name1 | grep -w " $atom " | head -n1 >> temp

    done

    ########################################################

    lines=`wc -l temp1 | awk '{ print $1 }'`
 
    if [ $lines -gt 0 ]; then
 
      for i in $(eval echo "{1..$lines}")
      do
        head -n $i temp1 | tail -n1 > temp2
        line1=`head -n1 temp2 | awk '{ print $0 }'`
        atomt=`head -n1 temp2 | awk '{ print $2 }'` 
        carga=`head -n1 temp2 | awk '{ print $7 }'`
        sed -i "s/$carga/ 0.0000/" temp2
        sed -i "s/  $atomt/${atomt}CV/" temp2
        line2=`head -n1 temp2 | awk '{ print $0 }'`
        sed -i "s/$line1/$line2/" $name1
      done
    fi
    sed -i "s/ XXCV / CCV\* /g" $name1
    sed -i "s/ YYCV / NCV\* /g" $name1
    sed -i "s/ XX / C\* /g" $name1
    sed -i "s/ YY / N\* /g" $name1

    sed -i "s/ C3RCV /  C3RC /g" $name1
    sed -i "s/ C2RCV /  C2RC /g" $name1

    ########################################################

#    lines=`wc -l temppp | awk '{ print $1 }'`
#
#    if [ $lines -gt 0 ]; then
#
#      for i in $(eval echo "{1..$lines}")
#      do
#        head -n $i temppp | tail -n1 > temp2
#        line1=`head -n1 temp2 | awk '{ print $0 }'`
#        atomt=`head -n1 temp2 | awk '{ print $2 }'`
#        sed -i "s/  $atomt/${atomt}CV/" temp2
#        line2=`head -n1 temp2 | awk '{ print $0 }'`
#        sed -i "s/$line1/$line2/" $name2
#      done
#    fi
  
#    grep -w " RET " $name2 | grep -w " CE " > temp
#    grep -w " RET " $name2 | grep -w " HE. " >> temp
#    grep -w " RET " $name2 | grep -w " NZ " >> temp
#    grep -w " RET " $name2 | grep -w " HZ1 " >> temp
#    grep -w " C3R \| C2R \| HR " $name2 >> temp

    lines=`wc -l temp | awk '{ print $1 }'`

    if [ $lines -gt 0 ]; then

      for i in $(eval echo "{1..$lines}")
      do
        head -n $i temp | tail -n1 > temp2
        line1=`head -n1 temp2 | awk '{ print $0 }'`
        atomt=`head -n1 temp2 | awk '{ print $2 }'`
        carga=`head -n1 temp2 | awk '{ print $7 }'`
        sed -i "s/$carga/ 0.0000/" temp2
        sed -i "s/  $atomt/${atomt}CV/" temp2
        line2=`head -n1 temp2 | awk '{ print $0 }'`
        sed -i "s/$line1/$line2/" $name2
      done
    fi

#    if [ $lines -gt 0 ]; then
#
#      for i in $(eval echo "{1..$lines}")
#      do
#        head -n $i temp | tail -n1 > temp2
#        line1=`head -n1 temp2 | awk '{ print $0 }'`
#        atomt=`head -n1 temp2 | awk '{ print $2 }'`
#        sed -i "s/  $atomt/${atomt}CV/" temp2
#        line2=`head -n1 temp2 | awk '{ print $0 }'`
#        sed -i "s/$line1/$line2/" $name3
#      done
#    fi

    sed -i "s/ XXCV / CCV\* /g" $name2
    sed -i "s/ YYCV / NCV\* /g" $name2
    sed -i "s/ XX / C\* /g" $name2
    sed -i "s/ YY / N\* /g" $name2
    sed -i "s/ C3RCV /  C3RC /g" $name2
    sed -i "s/ C2RCV /  C2RC /g" $name2


#    sed -i "s/ XX / C\* /g" $name3
#    sed -i "s/ YY / N\* /g" $name3

#    sed -i "s/ C3RCV /  C3RC /g" $name3
#    sed -i "s/ C2RCV /  C2RC /g" $name3

    rm temp temp1 temp2

  done

  fil=`ls -1 ${Project}_Protein_chain_?.itp | wc -l`

  if [ $fil -eq 1 ]; then
    cp ${Project}_Protein_chain_A_1.itp Rerun_1/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_A_2.itp Rerun_2/${Project}_Protein_chain_A.itp
#    cp ${Project}_Protein_chain_A_3.itp Rerun_3/${Project}_Protein_chain_A.itp
    cp ../../Dynamic/output/${Project}_Ion.itp Rerun_1
    cp ../../Dynamic/output/${Project}_Ion.itp Rerun_2
#    cp ../../Dynamic/output/${Project}_Ion.itp Rerun_3
  fi
  if [ $fil -eq 2 ]; then
    cp ${Project}_Protein_chain_A_1.itp Rerun_1/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_A_2.itp Rerun_2/${Project}_Protein_chain_A.itp
#    cp ${Project}_Protein_chain_A_3.itp Rerun_3/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_B_1.itp Rerun_1/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_B_2.itp Rerun_2/${Project}_Protein_chain_B.itp
#    cp ${Project}_Protein_chain_B_3.itp Rerun_3/${Project}_Protein_chain_B.itp
    cp ../../Dynamic/output/${Project}_Ion_chain_B2.itp Rerun_1
    cp ../../Dynamic/output/${Project}_Ion_chain_B2.itp Rerun_2
#    cp ../../Dynamic/output/${Project}_Ion_chain_B2.itp Rerun_3  
  fi
  if [ $fil -eq 3 ]; then
    cp ${Project}_Protein_chain_A_1.itp Rerun_1/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_A_2.itp Rerun_2/${Project}_Protein_chain_A.itp
#    cp ${Project}_Protein_chain_A_3.itp Rerun_3/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_B_1.itp Rerun_1/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_B_2.itp Rerun_2/${Project}_Protein_chain_B.itp
#    cp ${Project}_Protein_chain_B_3.itp Rerun_3/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_C_1.itp Rerun_1/${Project}_Protein_chain_C.itp
    cp ${Project}_Protein_chain_C_2.itp Rerun_2/${Project}_Protein_chain_C.itp
#    cp ${Project}_Protein_chain_C_3.itp Rerun_3/${Project}_Protein_chain_C.itp
    cp ../../Dynamic/output/${Project}_Ion_chain_C2.itp Rerun_1
    cp ../../Dynamic/output/${Project}_Ion_chain_C2.itp Rerun_2
#    cp ../../Dynamic/output/${Project}_Ion_chain_C2.itp Rerun_3    
  fi
  if [ $fil -eq 4 ]; then
    cp ${Project}_Protein_chain_A_1.itp Rerun_1/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_A_2.itp Rerun_2/${Project}_Protein_chain_A.itp
#    cp ${Project}_Protein_chain_A_3.itp Rerun_3/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_B_1.itp Rerun_1/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_B_2.itp Rerun_2/${Project}_Protein_chain_B.itp
#    cp ${Project}_Protein_chain_B_3.itp Rerun_3/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_C_1.itp Rerun_1/${Project}_Protein_chain_C.itp
    cp ${Project}_Protein_chain_C_2.itp Rerun_2/${Project}_Protein_chain_C.itp
#    cp ${Project}_Protein_chain_C_3.itp Rerun_3/${Project}_Protein_chain_C.itp
    cp ${Project}_Protein_chain_D_1.itp Rerun_1/${Project}_Protein_chain_D.itp
    cp ${Project}_Protein_chain_D_2.itp Rerun_2/${Project}_Protein_chain_D.itp
#    cp ${Project}_Protein_chain_D_3.itp Rerun_3/${Project}_Protein_chain_D.itp
    cp ../../Dynamic/output/${Project}_Ion_chain_D2.itp Rerun_1
    cp ../../Dynamic/output/${Project}_Ion_chain_D2.itp Rerun_2
#    cp ../../Dynamic/output/${Project}_Ion_chain_D2.itp Rerun_3
  fi
  if [ $fil -eq 5 ]; then
    cp ${Project}_Protein_chain_A_1.itp Rerun_1/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_A_2.itp Rerun_2/${Project}_Protein_chain_A.itp
#    cp ${Project}_Protein_chain_A_3.itp Rerun_3/${Project}_Protein_chain_A.itp
    cp ${Project}_Protein_chain_B_1.itp Rerun_1/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_B_2.itp Rerun_2/${Project}_Protein_chain_B.itp
#    cp ${Project}_Protein_chain_B_3.itp Rerun_3/${Project}_Protein_chain_B.itp
    cp ${Project}_Protein_chain_C_1.itp Rerun_1/${Project}_Protein_chain_C.itp
    cp ${Project}_Protein_chain_C_2.itp Rerun_2/${Project}_Protein_chain_C.itp
#    cp ${Project}_Protein_chain_C_3.itp Rerun_3/${Project}_Protein_chain_C.itp
    cp ${Project}_Protein_chain_D_1.itp Rerun_1/${Project}_Protein_chain_D.itp
    cp ${Project}_Protein_chain_D_2.itp Rerun_2/${Project}_Protein_chain_D.itp
#    cp ${Project}_Protein_chain_D_3.itp Rerun_3/${Project}_Protein_chain_D.itp
    cp ${Project}_Protein_chain_E_1.itp Rerun_1/${Project}_Protein_chain_E.itp
    cp ${Project}_Protein_chain_E_2.itp Rerun_2/${Project}_Protein_chain_E.itp
#    cp ${Project}_Protein_chain_E_3.itp Rerun_3/${Project}_Protein_chain_E.itp
    cp ../../Dynamic/output/${Project}_Ion_chain_E2.itp Rerun_1
    cp ../../Dynamic/output/${Project}_Ion_chain_E2.itp Rerun_2
#    cp ../../Dynamic/output/${Project}_Ion_chain_E2.itp Rerun_3
  fi

  cd ..

fi

cd Total_CASSCF
cp -r ../../Dynamic/amber94.ff Rerun_1
cp -r ../../Dynamic/amber94.ff Rerun_2
#cp -r ../../Dynamic/amber94.ff Rerun_3

cp ../../Dynamic/output/$Project.gro Rerun_1
#cp ../../Dynamic/output/${Project}_Ion.itp Rerun_1
cp ../../Dynamic/output/$Project.top Rerun_1
cp ../../Dynamic/output/$Project.ndx Rerun_1
cp ../../Dynamic/output/$Project.trr Rerun_1
cp ../../Dynamic/output/posre* Rerun_1
cp ../../Dynamic/output/residuetypes.dat Rerun_1
#cp ../../Dynamic/output/traj.xtc Rerun_1
cp ../../Dynamic/output/dynamic.mdp Rerun_1

cp ../../Dynamic/output/$Project.gro Rerun_2
#cp ../../Dynamic/output/${Project}_Ion.itp Rerun_2
cp ../../Dynamic/output/$Project.top Rerun_2
cp ../../Dynamic/output/$Project.ndx Rerun_2
cp ../../Dynamic/output/$Project.trr Rerun_2
cp ../../Dynamic/output/posre* Rerun_2
cp ../../Dynamic/output/residuetypes.dat Rerun_2
#cp ../../Dynamic/output/traj.xtc Rerun_2
cp ../../Dynamic/output/dynamic.mdp Rerun_2

#cp ../../Dynamic/output/$Project.gro Rerun_3
##cp ../../Dynamic/output/${Project}_Ion.itp Rerun_3
#cp ../../Dynamic/output/$Project.top Rerun_3
#cp ../../Dynamic/output/$Project.ndx Rerun_3
#cp ../../Dynamic/output/$Project.trr Rerun_3
#cp ../../Dynamic/output/posre* Rerun_3
#cp ../../Dynamic/output/residuetypes.dat Rerun_3
##cp ../../Dynamic/output/traj.xtc Rerun_3
#cp ../../Dynamic/output/dynamic.mdp Rerun_3


cp -f  $templatedir/ASEC/ffbonded_mod.itp Rerun_1/amber94.ff/ffbonded.itp
cp -f  $templatedir/ASEC/ffnonbonded_mod.itp Rerun_1/amber94.ff/ffnonbonded.itp
cp -f  $templatedir/ASEC/ffbonded_mod.itp Rerun_2/amber94.ff/ffbonded.itp
cp -f  $templatedir/ASEC/ffnonbonded_mod.itp Rerun_2/amber94.ff/ffnonbonded.itp
#cp -f  $templatedir/ASEC/ffbonded_mod.itp Rerun_3/amber94.ff/ffbonded.itp
#cp -f  $templatedir/ASEC/ffnonbonded_mod_3.itp Rerun_3/amber94.ff/ffnonbonded.itp

sol=0
sol=`grep -c " SOL " list_gro.dat`

if [ $sol -gt 0 ]; then
#   cp -f  $templatedir/ASEC/tip3p_1.itp Rerun_1/amber94.ff/tip3p.itp
   cp -f  $templatedir/ASEC/tip3p_2.itp Rerun_2/amber94.ff/tip3p.itp
#   cp -f  $templatedir/ASEC/tip3p_3.itp Rerun_3/amber94.ff/tip3p.itp

   grep " SOL " list_gro.dat > temp

   for i in $(eval echo "{1..$sol}")
   do
     solnum=`head -n $i temp | tail -n1 | awk '{ print $1 }'`
     cd Rerun_2
     solline=`grep "SOL " $Project.gro | grep "$solnum " | awk '{ print $0 }'`
     grep "SOL " $Project.gro | grep "$solnum " > temp2
     sed -i "s/SOL/SCV/" temp2
     solline2=`head -n1 temp2 | awk '{ print $0 }'`
     sed -i "s/$solline/$solline2/" $Project.gro
     rm temp2
     cd ..
   done
   rm temp
   cd Rerun_2
   topnum=`grep "SOL          " $Project.top | awk '{ print $2 }'`
   topnew1=$(($topnum-$sol/3))
   topnew2=$(($sol/3))
   sed -i "/SOL          /d" $Project.top
   echo "SOL                $topnew1" >> $Project.top
   echo "SCV                $topnew2" >> $Project.top
#   cp $Project.gro ../Rerun_2
#   cp $Project.top ../Rerun_2
#   cp $Project.gro ../Rerun_3
#   cp $Project.top ../Rerun_3
   cd ..
fi

cd Rerun_1
$gropath/gmx grompp -f dynamic.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr > mdlog
$gropath/gmx mdrun -rerun $Project.trr -s $Project.tpr
cd ..
cd Rerun_2
$gropath/gmx grompp -f dynamic.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr > mdlog
$gropath/gmx mdrun -rerun $Project.trr -s $Project.tpr
#cd ..
#cd Rerun_3
#$gropath/gmx grompp -f dynamic.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr > mdlog
#$gropath/gmx mdrun -rerun $Project.trr -s $Project.tpr

cd ..
cp Rerun_1/md.log md_1.log
cp Rerun_2/md.log md_2.log
#cp Rerun_3/md.log md_3.log
#cp ../../Dynamic/output/md.log .

cp ../${Project}_6-31G_Opt/${Project}_6-31G_Opt.out molcas.out 

cp $templatedir/ASEC/Energies_CAV_CASSCF.f .
#cp ../Energies_CAV.f .
gfortran Energies_CAV_CASSCF.f -o Energies_CAV_CASSCF.x
./Energies_CAV_CASSCF.x
rm Energies_CAV_CASSCF.x Energies_CAV_CASSCF.f

cd ../

###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

