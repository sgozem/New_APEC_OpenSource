#!/bin/bash
#
# Gromacs (all atoms) --> Tinker PDB conversion
#
# The PDB to use must be the 1st argument of the script
# If it is not found, the script terminates, otherwise some info is given to the user
#
if [ -z $1 ]; then
   echo ""
   echo " No PDB input file specified."
   echo " You have to type ./pdb-format.sh followed by the PDB file name"
   echo " Terminating..." 
   echo ""
   exit 0
fi
pdbfile=$1
if [ ! -f $pdbfile ]; then
   echo ""
   echo " PDB file not existing. Check if all files are in the same folder."
   echo " Terminating..."
   echo ""
   exit 0
else
   echo ""
   echo " $pdbfile will be modified and saved to final-tk.pdb"
   echo " final-tk.pdb is ready for Tinker pdbxyz conversion"
   echo ""
   echo " Now processing..."
   echo ""
fi
#
# Backing up $pdbfile
#
cp $pdbfile bak_$pdbfile
#
# Generation of replacements vim script
#
echo ':%s/1HH3/ H1 /' > replacements.vim
echo ':%s/2HH3/ H2 /' >> replacements.vim
echo ':%s/3HH3/ H3 /' >> replacements.vim
echo ':%s/OC1/O  /' >> replacements.vim
echo ':%s/OC2/O2 /' >> replacements.vim
echo ':g/SOL/s/ATOM  /HETATM/' >> replacements.vim
echo ':%s/SOL/HOH/' >> replacements.vim
echo ':g/HOH/s/OW/O /' >> replacements.vim
echo ':g/HOH/s/HW1/H1 /' >> replacements.vim
echo ':g/HOH/s/HW2/H2 /' >> replacements.vim
echo ':x' >> replacements.vim
vi -es $pdbfile < replacements.vim
rm replacements.vim

cp $pdbfile final-tk.pdb 

#
# Changing some hydrogen atoms labels
#
sed -i "s/1HD1/HD11/;s/2HD1/HD12/;s/3HD1/HD13/" final-tk.pdb
sed -i "s/1HD2/HD21/;s/2HD2/HD22/;s/3HD2/HD23/" final-tk.pdb
sed -i "s/1HG1/HG11/;s/2HG1/HG12/;s/3HG1/HG13/" final-tk.pdb
sed -i "s/1HG2/HG21/;s/2HG2/HG22/;s/3HG2/HG23/" final-tk.pdb
#
# Ions need to be HETATM
#
sed -i "/CL/s/ATOM  /HETATM/; /NA/s/ATOM  /HETATM/; /ACI/s/ATOM  /HETATM/" final-tk.pdb
#
# Generating the vim script for chain label replacement if multiple chains are detected
#
multchain=`grep "MultChain" ../Infos.dat | awk '{ print $2 }'`
if [[ $multchain != "YES" && $multchain != "NO" ]]; then
   echo " Cannot read multiple chains information from Infos.dat"
   echo " Are there multiple chains in your protein? (y/n)"
   read answer
   echo ""
   if [[ $answer == "y" ]]; then
      multchain=YES
   else
      multchain=NO
   fi
fi
if [[ $multchain == "YES" ]]; then
   IFS=" "
   lastresa=( $( grep "LastRes" ../Infos.dat ) )
   diffchain=( $( grep "DiffChain" ../Infos.dat ) )
   ngap=${#lastresa[@]}-1
   lastgap=$(($ngap+1))
   lastresa[$lastgap]=1000
   i=1
   for lette in {A..E}; do
       catene[$i]=$lette
       i=$(($i+1))
   done
   IFS=""
   for ((i=1;i<$(($lastgap));i=$(($i+1)))); do
       j=$(($i+1))
       awk '{ if ( $6 > '"${lastresa[$i]}"' ) sub(" '"${catene[$i]}"' "," '"${catene[$j]}"' "); print $0 }' final-tk.pdb > new2.pdb
       mv new2.pdb final-tk.pdb
   done
fi

#
# Temporary files are removed. Comment these 2 lines for debug
#
cp bak_$pdbfile $pdbfile 
rm bak_$pdbfile
#
# Normal termination
# 
echo " Normal termination - final-tk.pdb was created successfully"
echo " Have a nice day!"
echo ""

