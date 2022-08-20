#!/bin/bash
#
#
# Changing polar hydrogens labels to match Gromacs standard after Dowser execution
# Called by rundowser.sh, with pdb filename as the only argument
#
Project=$1
#
# Isoleucine C delta, Threonine H gamma, Arginine hydrogens
#
sed "/ILE/s/CD1/CD /" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/THR/s/HG /HG1/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/ARG/s/2HH2/HH22/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/ARG/s/2HH1/HH12/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/ARG/s/1HH1/HH11/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/ARG/s/1HH2/HH21/" $Project.pdb > temporal
mv temporal $Project.pdb
#
# Asparagine, glutamine
#
sed "s/ HD1 ASN/HD21 ASN/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "s/ HD2 ASN/HD22 ASN/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "s/ HE1 GLN/HE21 GLN/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "s/ HE2 GLN/HE22 GLN/" $Project.pdb > temporal
mv temporal $Project.pdb
#
# Neutral Asp H delta. Why not GLH?
#
sed "/ASH/s/HD /HD2/" $Project.pdb > temporal
mv temporal $Project.pdb
#
# Messages to the user
#
echo " Hydrogen atoms labels changed successfully, now modifying water and histidines"
echo ""
#
# Changing labels of water molecules
#
sed "/HOH/s/H1 /HW1/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "/HOH/s/H2 /HW2/" $Project.pdb > temporal
mv temporal $Project.pdb
sed "s/HOH W/HOH A/" $Project.pdb > temporal
mv temporal $Project.pdb

