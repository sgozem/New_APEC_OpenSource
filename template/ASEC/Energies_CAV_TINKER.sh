#!/bin/bash
#
# Retrieving information from Infos.dat
#

Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
mode=`grep "Mode " ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" ../Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" ../Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" ../Infos.dat | awk '{ print $2 }'`

config=1000

##
##  Here will be generated the all the tinker files from 
##  the selected configurations of the molecular dynamic

mkdir Total_TINKER
cd Total_TINKER
mkdir Tinker_calc
cd Tinker_calc
##########################
#generating the keyfile using the list_tk.dat
#########################

cp ../../${Project}_6-31G_Opt/${Project}_6-31G_Opt.key Tinker_cav.key
sed -i "/ACTIVE/d" Tinker_cav.key

cavnum=`grep "atoms" ../../../MD_analysis/list_tk.dat | awk '{print $4}'`

for i in $(eval echo "{1..$cavnum}")
do
  cavatom=`head -n $(($i+1)) ../../../MD_analysis/list_tk.dat | tail -n1 | awk '{print $1}'`
  echo "ACTIVE $cavatom" >> Tinker_cav.key
done

cont=`grep -w -c "MM\|QM\|LA" Tinker_cav.key`

for i in $(eval echo "{1..$cont}")
do
  temp2=`grep -w -m $i "MM\|QM\|LA" Tinker_cav.key | tail -n1 | awk '{print $2}'`
  temp3=`grep -w -m $i "MM\|QM\|LA" Tinker_cav.key | tail -n1 | awk '{print $3}'`
  echo "ACTIVE $temp2 $temp3" >> Tinker_cav.key
done
cp Tinker_cav.key $Project-tk.key

sed -i "/LA /d" $Project-tk.key
lanum=`grep "QMMM " $Project-tk.key | awk '{print $2}'`
sed -i "s|QMMM $lanum|QMMM $(($lanum-1))|" $Project-tk.key

######################################################################

#
# This section is for inserting the optimizaed cromophore structure of the previous 
# Step inside the new -tk.xyz in order to avoid the rounding arising from the tk to gro
# conversion and the new position of the Link atom.
###   (Taken from Molcami_CASSCF_oneIter.sh) ###

cp $templatedir/ASEC/Ins_tk_tk_Molcami.f .
cp ../../../${Project}_6-31G_Opt.Final_last.xyz last_tk.xyz

cont2=`grep -w -c "MM\|QM" Tinker_cav.key`

sed -i "s|intervalos|${cont2}|" Ins_tk_tk_Molcami.f
sed -i "s|numero|$(($numatoms+1))|" Ins_tk_tk_Molcami.f

for i in $(eval echo "{$cont2..1}")
do
   temp4=`grep -w -m $i "MM\|QM" Tinker_cav.key | tail -n1 | awk '{print $3}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,2)=$temp4" Ins_tk_tk_Molcami.f
   temp5=`grep -w -m $i "MM\|QM" Tinker_cav.key | tail -n1 | awk '{print ((-1*$2))}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,1)=$temp5" Ins_tk_tk_Molcami.f
done

gfortran Ins_tk_tk_Molcami.f -o Ins_tk_tk_Molcami.x

#cp ../../../Dynamic/output/$Project.tpr .
#cp ../../../Dynamic/output/traj.xtc .
$gropath/gmx trjconv -s ../../../Dynamic/output/$Project.tpr -f ../../../Dynamic/output/traj.xtc -b 301 -skip 1 -o Selected.gro

#cp ../../../MD_analysis/Selected_40.gro Selected.gro
cp ../../${Project}_6-31G_Opt/$prm.prm .
cp $templatedir/ASEC/pdb-format-new_mod.sh .

##   (This section is from MD_2_QM.sh)

#=1
for j in $(eval echo "{1..$config}")
do
   head -n $((($numatoms+3)*$j)) Selected.gro | tail -n $(($numatoms+3)) > final-${Project}.gro

   # editconf convert it to PDB and pdb-format-new fixes the format to
   # allow Tinker reading
         
#   cp $templatedir/ASEC/pdb-format-new_mod.sh .
   $gropath/gmx editconf -f final-${Project}.gro -o final-$Project.pdb -label A
   echo "Configuration $j"
   ./pdb-format-new_mod.sh final-$Project.pdb
   # pdbxyz conversion
   #            
   mv final-tk.pdb $Project-tk.pdb
   $tinkerdir/pdbxyz $Project-tk.pdb << EOF
ALL
$prm
EOF

#
# xyzedit generates the Molcas input file and adds the link atom
#

   $tinkerdir/xyzedit $Project-tk.xyz <<EOF
20
1
EOF

   rm final-${Project}.gro final-$Project.pdb $Project-tk.pdb $Project-tk.seq $Project-tk.xyz
   mv $Project-tk.xyz_2 new_tk.xyz
   rm $Project-tk.input

   ./Ins_tk_tk_Molcami.x
   mv Insert-tk.xyz Tinker_cav_$j.xyz
   rm new_tk.xyz

   cp $templatedir/ASEC/template_CAV Tinker_cav_$j.input
   cp Tinker_cav.key Tinker_cav_$j.key
done

rm $Project-tk.key last_tk.xyz
#cp $templatedir/ASEC/molcas-job_CAV.sh .

NOME=${Project}_$Step
#echo ""
#echo ""
#echo "Enter the Name of the Script"
#echo ""
#echo ""
#read NOME
for j in {1..10}
do
  cp $templatedir/ASEC/molcas-job_CAV_$j.sh .
  sed -i "s|NOMEPROGETTO|${NOME}_$j|" molcas-job_CAV_$j.sh
  qsub molcas-job_CAV_$j.sh
done


cp $templatedir/ASEC/Total_TINKER.sh ../../

echo ""
echo " After finishing the energy calculations run Total_TINKER.sh"
echo ""

