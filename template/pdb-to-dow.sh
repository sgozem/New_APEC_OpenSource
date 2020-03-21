#!/bin/bash
#
# PDB Tinker (all atoms) --> Dowser
#
# PDB file must be in the same directory, and its complete name must be given as an argument
# The directory of execution must include the file labelret, with proper retinal labels
# Called by rundowser.sh
#
# Retrieving retinal chain label and its lysine number from Infos.dat
#
pdb=$1
#retchain=`grep "RetChain" ../../Infos.dat | awk '{ print $2 }'`
#Lys=`grep "LysNum" ../../Infos.dat | awk '{ print $2 }'`
Project=`grep "Project" ../../Infos.dat | awk '{ print $2 }'`
#
# Checking for PDB file existence
#
if [ -z $pdb ]; then
   echo " File PDB not found! Aborting..."
   echo ""
   echo "pdb-to-dow.sh 1 PDBnotFound" >> ../../arm.err
   exit 0 
fi
# 
# Checking for labelret file existence
#
if [ -z labelret ]; then
   echo " File labelret not found!" 
   echo " Retinal labels file is needed! Aborting..."
   echo ""
   echo "pdb-to-dow.sh 2 NoLabelret" >> ../../arm.err
   exit 0
fi
#lyslabel=`grep " $retchain $Lys " $pdb | head -n 1 | awk '{ print $4 }'`
#
# Messages to the user
#
echo ""
echo " Converting $pdb Tinker PDB into Dowser format..."
echo ""
#
# Important residue numbers
#
#PLys=`grep -B1 "$lyslabel $retchain $Lys" $pdb | head -n 1 | awk '{ print $4 }'`
#
# Retinal residue number is given by the residue label of the first RET atom
#
#Ret=`grep 'RET' $pdb | head -n 1 | awk '{ print $6 }'`
# 
# Finding correct water labels, then 
# saving water molecules and chloride ions
#
grep 'HOH' $pdb > HOH
grep 'CL ' $pdb > CL
grep 'NA ' $pdb > NA
grep 'ACI' $pdb > ACI
#
# Saving retinal bound lysine in file LYS
# LYS is replaced by RET because in Gromacs the RET residue includes the lysine too
#
#grep "$lyslabel $retchain $Lys" $pdb  > LYS
#sed -i "1,$ s/$lyslabel/RET/" LYS
echo -e 'END' >> END
#
# The PDB file copy dowser.pdb is modified, so that the original PDB is preserved
#
cp $pdb dowser.pdb
#
# Removing garbage from PDB: CONECT records and MASTER. Also, non-retinal HETATM labels are deleted
#
sed '/CONECT/d;/MASTER/d' dowser.pdb > temporal
mv temporal dowser.pdb
sed '/RET/s/HETATM/ATOM  /;/HETATM/d' dowser.pdb > temporal
mv temporal dowser.pdb
#
# togliEND is the VI script to remove LYS atoms from dowser.pdb
#
#echo -e ':$' >> togliEND
#echo -e ':.d' >> togliEND
#echo -e ":g/$lyslabel $retchain $Lys/d" >> togliEND
#echo -e ':x' >> togliEND
#vi -es dowser.pdb < togliEND
#
# Lysine is appended to the end of dowser.pdb, right after RET, followed by Cl and water
# The residue number of the new RET residue is changed back to the lysine one
# to keep the residue number order
#
#cat LYS >> dowser.pdb 
cat HOH >> dowser.pdb
cat END >> dowser.pdb

#sed -i s/"RET $retchain $Ret"/"RET $retchain $Lys"/ dowser.pdb
#rm LYS togliEND
#
# Replacing the retinal labels with correct ones,
# and saving the new retinal in newret
# possible issue: column limits are hard-coded
#
#grep 'RET' dowser.pdb > onlyret
#echo -e ':$' >> change
#echo -e ':1,20 co .' >> change
#echo -e ':1,20d' >> change
#echo -e ':x' >> change 
#vi -es  onlyret < change
# 
# Checking if labelret has the same number of lines as testaret and codaret
#
#cut -c1-11 onlyret > testaret
#cut -c18-54 onlyret > codaret
#ntestal=`wc -l testaret | awk '{ print $1 }'`
#ncodal=`wc -l codaret | awk '{ print $1 }'`
#nlabel=`wc -l labelret | awk '{ print $1 }'`
#if [[ $ntestal -eq $nlabel && $ncodal -eq $nlabel ]]; then
#   paste -d" " testaret labelret codaret > newret
#   rm onlyret testaret codaret change
#   sed -i '/RET/d' dowser.pdb
#else
#   echo " Something wrong with retinal labels and/or labelret file! Aborting..."
#   echo ""
#   echo "pdb-to-dow.sh 3 RetLabelMismatch" >> ../../arm.err
#   exit 0
#fi
#
# baboomba is the VI script to place RET from the bottom of dowser.pdb back to where its lysine was
#
#PNum=$(($Lys-1))
#echo -e ":?$PLys $retchain $PNum" >> baboomba
#echo -e ':r newret' >> baboomba
#echo -e ':x' >> baboomba
#vi -es dowser.pdb  < baboomba
#rm baboomba newret
echo "pdb-to-dow.sh 0 OK" >> ../../arm.err
