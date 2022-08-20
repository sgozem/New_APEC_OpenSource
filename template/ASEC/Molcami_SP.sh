#!/bin/bash
#
# Retrieving the needed information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
uphess=`grep "Update_hessian" Infos.dat | awk '{ print $2 }'`

#
# Copying and running the script for Tinker files preparation
#

# Creating the directory calculations and copying all the needed files
# If calculations already exists, the script exits with a message
#
echo ""
echo " Using $Project.xyz and $prm.prm"
echo ""
if [ -d calculations ]; then
   echo " Folder calculations already exists. Please check it and remove if necessary"
   echo " Terminating..."
   echo ""
   exit 0
fi

cp $templatedir/ASEC/keymaker_mod.sh .
./keymaker_mod.sh $Project-tk $prm.prm

mkdir calculations
newdir=${Project}_6-31G
mkdir calculations/${newdir}

mv $Project-tk.key ${newdir}.key
sed -i "/ACTIVE/d" ${newdir}.key

#if [[ $Step -eq 0 ]]; then
#   cat ${newdir}.key Templates/list_tk > te
#   mv te ${newdir}.key
#fi

cont=`grep -w -c "MM\|QM\|LA" ${newdir}.key`

for i in $(eval echo "{1..$cont}")
do
temp2=`grep -w -m $i "MM\|QM\|LA" ${newdir}.key | tail -n1 | awk '{print $2}'`
temp3=`grep -w -m $i "MM\|QM\|LA" ${newdir}.key | tail -n1 | awk '{print $3}'`
echo "ACTIVE $temp2 $temp3" >> ${newdir}.key
done

#
# This section is for inserting the optimizaed cromophore structure of the previous 
# Step into the new -tk.xyz in order to avoid the rounding arising from the tk to gro
# conversion and the new position of the Link atom.
#
answer=b
if [[ $Step -eq 0 ]]; then
   while [[ $answer != y && $answer != n ]]; do
      echo " If you are optimizing a TS may be a good idea to insert a previously"
      echo " optimized chromophore structure into the -tk.xyz in order to avoid"
      echo " the rounding arising from the tinker to gromacs conversion."
      echo " If this is a minimum optimization you can skip this message."
      echo ""
      echo " Do you want to insert it? (y/n)"
      echo ""
      read answer
      if [[ $answer == "y" ]]; then
         echo " Enter the full path of the xyz file, including the name"
         read path
         cp $path ${newdir}.Final_last.xyz
      fi
   done
fi

if [[ $answer == "b" || $answer == "y" ]]; then
#if [ -f ${newdir}.Final_last.xyz ]; then
   cp $templatedir/ASEC/Ins_tk_tk_Molcami.f .
   cp $Project-tk.xyz new_tk.xyz
   cp ${newdir}.Final_last.xyz last_tk.xyz

   cont2=`grep -w -c "MM\|QM" ${newdir}.key`

   sed -i "s|intervalos|${cont2}|" Ins_tk_tk_Molcami.f
   sed -i "s|numero|$(($numatoms+1))|" Ins_tk_tk_Molcami.f

   for i in $(eval echo "{$cont2..1}")
   do
   temp4=`grep -w -m $i "MM\|QM" ${newdir}.key | tail -n1 | awk '{print $3}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,2)=$temp4" Ins_tk_tk_Molcami.f
   temp5=`grep -w -m $i "MM\|QM" ${newdir}.key | tail -n1 | awk '{print ((-1*$2))}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,1)=$temp5" Ins_tk_tk_Molcami.f
   done

   gfortran Ins_tk_tk_Molcami.f -o Ins_tk_tk_Molcami.x
   ./Ins_tk_tk_Molcami.x
   mv Insert-tk.xyz $Project-tk.xyz
   rm new_tk.xyz last_tk.xyz Ins_tk_tk_Molcami.f Ins_tk_tk_Molcami.x
fi

###################################################

##### cp new_$prm.prm calculations/$prm.prm

cp $Project-tk.xyz calculations/${newdir}/${newdir}.xyz
mv ${newdir}.key calculations/${newdir}/${newdir}.key
#cp $Project-tk.input calculations/${newdir}/${newdir}.input
#rm $Project-tk.input
#slurm
#cp $templatedir/molcas.slurm.sh calculations/${newdir}/molcas-job.sh
cp $templatedir/molcas-job.sh calculations/${newdir}/

if [[ $retstereo == "nAT" ]]; then
   cp $templatedir/templateOPTSCFneu calculations/${newdir}/templateSP
else
   cp $templatedir/ASEC/templateSP calculations/${newdir}/
fi

cd calculations
cd ${newdir}/

# Putting project name, input directory, time and memory in the submission script
# Since it is a SCF optimization, 1500 MB should be enough...
# And 30 hrs are enough as well
#

NOME=${newdir}
#echo ""
#echo ""
#echo "Enter the Name of the Script"
#echo ""
#echo ""
#read NOME
NOME=${Project}_6-31G
sed -i "s|NOMEPROGETTO|$NOME|" molcas-job.sh

no=$PWD
#sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|12000|" molcas-job.sh
sed -i "s|MEMORIA|10000|" molcas-job.sh
sed -i "s|hh:00:00|30:00:00|" molcas-job.sh
sed -i "s|ppn=4|ppn=4|" molcas-job.sh

# Replacing PARAMETRI with current prm filename templateOPTSCF
#
sed -i "s|PARAMETRI|${prm}|" templateSP

#slurm
#mv molcas-job.sh molcas.slurm.sh
sed -i "/#PBS -l mem=/a#PBS -A PAA0009" molcas-job.sh
sed -i "/export Project=/c\ export Project=${newdir}" molcas-job.sh

##echo "cp \$WorkDir/\$Project.JobIph \$InpDir/\$Project.JobIph_new" >> molcas-job.sh

mv templateSP ${newdir}.input 
sed -i "s|3-21G|6-31G*|" ${newdir}.input

cd ../../
cp $templatedir/ASEC/ASEC_SP.sh .
echo ""
echo "Run ASEC_SP.sh to generate the final coordinate file and submitt"
echo ""

