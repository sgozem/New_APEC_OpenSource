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
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`


  echo ""
  echo ""
  echo " Processing "
  echo ""
  echo ""

  cd MD_analysis
#  mkdir AutoCorrelation
  cd AutoCorrelation
  mkdir Rerun

  cp new_Dynamic/${Project}_Protein_chain_?.itp Rerun
#  cp ../../MD_analysis/list_gro.dat .
  
  cd Rerun
#######################
# In the Following the atom type of the Retinal atoms 
# will be changed to CV and chargees equal zero
#######################

  for k in ${Project}_Protein_chain_?.itp
  do

### ${k%.itp} get the basename of the files .itp, for instance:
### chain_A.itp become chain_A_1.itp and chain_A_2.itp
 
#    name1=${k%.itp}_1.itp
#    name2=${k%.itp}_2.itp
#    name3=${k%.itp}_3.itp

    name1=$k

#    cp $k $name1
#    cp $k $name2
#    cp $k $name3

    sed -i 's/ C\* / XX /g' $name1
    sed -i 's/ N\* / YY /g' $name1
#    sed -i 's/ C\* / XX /g' $name2
#    sed -i 's/ N\* / YY /g' $name2
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

    ########################################################

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
        sed -i "s/$line1/$line2/" $name1
      done
    fi
    sed -i "s/ XXCV / CCV\* /g" $name1
    sed -i "s/ YYCV / NCV\* /g" $name1
    sed -i "s/ XX / C\* /g" $name1
    sed -i "s/ YY / N\* /g" $name1

    sed -i "s/ C3RCV /  C3RC /g" $name1
    sed -i "s/ C2RCV /  C2RC /g" $name1

    if [ -f temp ]; then
       rm temp
    fi
    if [ -f temp2 ]; then
       rm temp2
    fi
    
  done

cp -r ../new_Dynamic/amber94.ff .

cp ../new_Dynamic/output/$Project.gro .
cp ../new_Dynamic/output/${Project}_Ion*.itp .
cp ../new_Dynamic/output/$Project.top .
cp ../new_Dynamic/output/$Project.ndx .
cp ../new_Dynamic/output/$Project.trr .
cp ../new_Dynamic/output/posre* .
cp ../new_Dynamic/output/residuetypes.dat .
#cp ../../Dynamic/output/traj.xtc Rerun_1
cp ../new_Dynamic/output/dynamic.mdp .
cp ../new_Dynamic/gromacs.sh .

cp -f $templatedir/ASEC/ffbonded_mod.itp amber94.ff/ffbonded.itp
cp -f $templatedir/ASEC/ffnonbonded_mod.itp amber94.ff/ffnonbonded.itp

$gropath/gmx grompp -f dynamic.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr > mdlog

sed -i "s|export inpdir=.*|export inpdir="$(pwd)"|g" gromacs.sh
sed -i "s|export outdir=.*|export outdir="$(pwd)"|g" gromacs.sh
sed -i "s|mdrun_mpi .*|mdrun_mpi -rerun "$Project".trr -s "$Project".tpr|" gromacs.sh

qsub gromacs.sh

echo ""
echo " Wait for rerun to end, then run Autocorrelation_3.sh "
echo ""

#cp md.log ../md_zero.log
#cp ../new_Dynamic/output/md.log ../
#cd ../

#cp $templatedir/ASEC/Autocorrelation.f .

#gfortran Autocorrelation.f -o Autocorrelation.x
#./Autocorrelation.x
#rm Autocorrelation.x

###############################################
###############################################

echo ""
echo ""
echo " Done  "
echo ""
echo ""

