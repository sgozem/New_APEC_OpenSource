#!/bin/bash
#
# Tinker (all atoms) --> Gromacs
#
# PDB file must be in the same directory, and its complete name must be given as an argument
# The directory of execution must include the file labelret, with proper retinal labels
#
# Retrieving retinal chain label and its lysine number from Infos.dat
# dowser
#
pdb=$1
dowser=$2
multchain=`grep "MultChain" ../Infos.dat | awk '{ print $2 }'`
#
# Checking for PDB file existence
#
if [ -z $pdb ]; then
   echo " File PDB not found! Aborting..."
   echo ""
   exit 0 
fi
#
# Messages to the user
#
echo ""
echo " Converting $pdb Tinker PDB into Gromacs format..."
echo ""
#
# From here, if dowser was used, most of the commands were executed before...
#
if [[ $dowser == "NO" ]]; then
# 
#  Finding correct water labels, then saving water molecules and chloride ions
#
   aqua=`grep "HOH" $pdb | head -n 1`
   if [[ -z $aqua ]]; then
      aqua='WAT'
   else
      aqua='HOH'
   fi
   grep $aqua $pdb > HOH
   grep 'CL ' $pdb > CLO
   grep 'NA ' $pdb > NAT
   echo -e 'END' >> HOH
fi
#
# The PDB file copy new.pdb is modified, so that the original PDB is preserved
#
cp $pdb new.pdb
#
# The following is required only when dowser is not used
#
if [[ $dowser == "NO" ]]; then
#
#  Removing garbage from PDB: CONECT records and MASTER. Also, non-retinal HETATM labels are deleted
#
   sed -i '/CONECT/d;/MASTER/d' new.pdb
   sed -i '/RET/s/HETATM/ATOM  /;/HETATM/d' new.pdb
#
#  togliEND is the VI script to remove water, Cl ions, and remove LYS atoms
#  from new.pdb
#
   echo -e ":g/$aqua/d" >> togliEND
   echo -e ':$' >> togliEND
   echo -e ':.d' >> togliEND
   echo -e ':x' >> togliEND
   vi -es new.pdb < togliEND
#
#  Lysine is appended to the end of new.pdb, right after RET, followed by Cl and water
#  The residue number of the new RET residue is changed back to the lysine one
#  to keep the residue number order
#
   cat HOH >> new.pdb
   sed -i s/"RET $retchain $Ret"/"RET $retchain $Lys"/ new.pdb
   rm togliEND HOH
fi
#
# Removing remarks
#
sed -i "/REMARK/d" new.pdb
#
# Gromacs requires different labeling for the last residue
# OT labels is required for Gromacs 4.5.5 if DOWSER is YES ??
# In case there are problems with other version of GROMACS, the following might need changes
#
if [[ $dowser == "YES" ]]; then
   oxyterm="OT"
#### YOELVIS, moved here from line 172 (see below) for consistency
# change labels from OT to OXT for GROMACS 5.0.6 (doesn't affect 4.5.5)
   sed -i "s/$oxyterm /OXT/" new.pdb
#########
else
   oxyterm="OXT"
fi
if [[ $multchain == "YES" ]]; then
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
   for ((i=1;i<$(($lastgap));i=$(($i+1)))); do
       j=$(($i+1))
       awk '{ if ( $6 > '"${lastresa[$i]}"' ) sub(" '"${catene[$i]}"' "," '"${catene[$j]}"' "); print $0 }' new.pdb > new2.pdb
       mv new2.pdb new.pdb
   done
# YOE moved up
# change labels from OT to OXT for GROMACS 5.0.6 (doesn't affect 4.5.5)
#   sed -i "s/$oxyterm /OXT/" new.pdb
else
   lastres=`grep $oxyterm new.pdb | head -n 1 | awk '{ print $6 }'`
   lastype=`grep $oxyterm new.pdb | head -n 1 | awk '{ print $4 }'`
   sed -i "s/$oxyterm/OC2/" new.pdb
   sed -i "s/O   $lastype A $lastres/OC1 $lastype A $lastres/" new.pdb
fi

#
# Changing the label for protonated histidines
#
sed -i "s/HIS/HIP/" new.pdb
#
# Chloride ions must be added only when dowser is not used
#
if [[ $dowser == "NO" ]]; then
   cat CLO >> new.pdb
   cat NAT >> new.pdb
   rm CLO NAT
fi
echo "pdb-to-gro.sh 0 OK" > ../arm.err

