#!/bin/bash
#
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
dowser=`grep "Dowser" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
chr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`
chargechr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`

echo ""
echo " Name of the project is ${Project}"
echo " ${prm}.prm file found, and i am using it..."
echo ""
#
# Creating the directory where hydrogen atoms MM minimization will be done
# and putting the required files
#
if [[ -d Minimize_${Project} ]]; then
   ./smooth_restart.sh Minimize_${Project} "Do you want to re-run Dowser + H minimization? (y/n)" 1
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir Minimize_${Project}
cp $Project.pdb Minimize_${Project}
cd Minimize_${Project}
cp -r $templatedir/$amber.ff .

cp $templatedir/residuetypes.dat .
cp $templatedir/standard-EM.mdp .
cp $templatedir/soglia .
cp $templatedir/pdb-to-gro.sh .

#
# Befor the energy minimization dowser will be used to deal with external/internal waters
#
if [[ $dowser == "YES" ]]; then
   mkdir ${Project}_dowser
   cp $templatedir/carbret_dow ${Project}_dowser/labelret
   cp $templatedir/pdb-to-dow.sh ${Project}_dowser/
   cp $templatedir/yesH-tk-to-gro.sh ${Project}_dowser/
   cp $templatedir/${prm}.prm ${Project}_dowser/
   cp ../$Project.pdb ${Project}_dowser
   cp $templatedir/PdbFormatter.py ${Project}_dowser
   cp $templatedir/rundowser.sh ${Project}_dowser/
   cd ${Project}_dowser/
   ./rundowser.sh $Project $tinkerdir $prm
   checkrundow=`grep rundowser ../../arm.err | awk '{ print $2 }'`
   if [[ $checkrundow -ne 0 ]]; then
      echo " Problem in rundowser.sh. Aborting..."
      echo ""
      echo " NewStep.sh 1 RunDowserProbl" >> ../../arm.err
      exit 0
   fi
   mv ../$Project.pdb ../$Project.pdb.old.1
   cp $Project.pdb ../
   cd ..
######
#  Dowser is adding HZ1 even though nAT is selected
#   if [[ $retstereo == "nAT" ]]; then
#      sed -i "/ HZ1 RET /d" $Project.pdb
#   fi
######
else
   cp $templatedir/carbret labelret
fi

#
# Converting the PDB into a format suitable for Gromacs by using pdb-to-gro.sh
# Output must be different if Dowser was used
#
./pdb-to-gro.sh $Project.pdb $dowser
if [[ -f ../arm.err ]]; then
   checkpdbdow=`grep 'pdb-to-gro' ../arm.err | awk '{ print $2 }'`
fi
if [[ $checkpdbdow -ne 0 ]]; then
   echo " An error occurred in pdb-to-gro.sh. I cannot go on..."
   echo ""
   echo "NewStep.sh 2 PDBtoGroProblem" >> ../arm.err
   exit 0
fi
#
# Backing up the starting PDB and renaming new.pdb (the output of pdb-to-gro.sh)
#
mv $Project.pdb $Project.pdb.old.2
mv new.pdb $Project.pdb
echo " $Project.pdb converted successfully! Now it will converted into $Project.gro, "
echo " the Gromacs file format"
echo ""

wat=`grep -c "OW  HOH" $Project.pdb`
../update_infos.sh "DOWSER_wat" $wat ../Infos.dat

#
# pdb2gmx is the Gromacs utility for generating gro files and topologies
#

$gropath/gmx pdb2gmx -f $Project.pdb -o $Project.gro -p $Project.top -ff $amber -water tip3p 2> grolog
checkgro=`grep 'Writing coordinate file...' grolog`
   if [[ -z $checkgro ]]; then
      echo " An error occurred during the execution of pdb2gmx. Please look into grolog file"
      echo " No further operation performed. Aborting..."
      echo ""
      exit 0
   else
      echo " new.gro and its topology were successfully generated"
      echo ""
      rm grolog
   fi

#
# The tampletes gro2tk and tk2gro will be created for future
# back and forwards conversions
#
echo " **********************************************************************"
echo ""
echo " The Gromacs - Tinker templates will be created"
echo ""
echo " **********************************************************************"
cd ..

mkdir Templates
cp Minimize_$Project/$Project.gro Templates
cp $templatedir/$prm.prm Templates
cd Templates
sed -i "s/HOH/SOL/g" $Project.gro

#
# editconf convert it to PDB and pdb-format-new fixes the format to
# allow Tinker reading
#
cp $templatedir/ASEC/pdb-format-new_mod.sh .
$gropath/gmx editconf -f ${Project}.gro -o $Project.pdb -label A
./pdb-format-new_mod.sh $Project.pdb

#
# pdbxyz conversion
#
mv final-tk.pdb $Project-tk.pdb

# If PRO is a terminal residue (N-terminal or residue 1) the extra hydrogen is labeled in 
# GROMACS as H2, being H1 and H2 the hydrogens bonded to the N. But in TINKER
# (specifically in the pdbxyz) these hydrogens are labeled as H2 and H3. So, it will be relabeled.
# This is also performed in MD_2_QMMM.sh
sed -i "s/ATOM      3  H2  PRO A   1 /ATOM      3  H3  PRO A   1 /" $Project-tk.pdb
sed -i "s/ATOM      2  H1  PRO A   1 /ATOM      2  H2  PRO A   1 /" $Project-tk.pdb

$tinkerdir/pdbxyz $Project-tk.pdb << EOF
ALL
$prm
EOF
echo " Please wait ..."

numatoms=`head -n2 $Project.gro | tail -n1 | awk '{ print $1 }'`

cp $templatedir/ASEC/Templates_gro_tk.f .
sed -i "s|numero|$numatoms|g" Templates_gro_tk.f

cp $Project.gro final_Config.gro
cp $Project-tk.xyz coordinates_tk.xyz
#   cp $templatedir/residuetypes.dat .
gfortran Templates_gro_tk.f -o Templates_gro_tk.x
./Templates_gro_tk.x
rm final_Config.gro
rm coordinates_tk.xyz
cp template_tk2gro ../
cp template_gro2tk ../

cp $templatedir/ASEC/Solvent_box.sh ../
../update_infos.sh "Next_script" "Solvent_box.sh" ../Infos.dat

echo ""
echo "**********************************************************"
echo " "
echo " Continue executing: Solvent_box.sh"
echo ""
echo "**********************************************************"
echo ""

