#!/bin/bash
#
#
# The first argument is a PDB file, the second is what to move:
# 1 - Hydrogen within 7 Angstroms
# 2 - Sidechains within 4 Angstroms
#
pdb=$1.gro
moving=$2
backbone=$3
risposta=$4
#
# Reading the distance within sidechains will be moved or if the cavity file is available
#
if [[ -f ../Infos.dat ]]; then
   radius=`grep "RadiusMD" ../Infos.dat | awk '{ print $2 }'`
   Hraggio=`grep "RadiusH" ../Infos.dat | awk '{ print $2 }'`
   cavity=`grep "CavityFile" ../Infos.dat | awk '{ print $2 }'`
else
   echo " No Infos.dat found! The radius for sidechain will be set to the default value of 4"
   echo " and that for hydrogens will be set to 7"
   echo ""
   radius=4
   Hraggio=7
fi
if [[ -z $radius ]]; then
   radius=4
fi
if [[ -z $Hraggio ]]; then
   Hraggio=7
fi
#
# If a cavity file is found, the array of residue IDs is built
#
if [[ $cavity == YES ]]; then
   if [[ -f ../cavity ]]; then
      echo "resid" > testa
      cat testa ../cavity > cavitysel
      residline=`tr '\n' ' ' < cavitysel`
      rm testa cavitysel
   else
      echo " WARNING! Cavity file should be present according to Infos.dat, but it cannot be found!"
      echo " Going on with default distance-based cavity selection:"
      echo " $radius for sidechains and $Hraggio for hydrogens"
      echo ""
      cavity=NO
   fi
fi

#
# Checking for PDB file existence
#
if [ -f $pdb ]; then
   echo " Creating the ndx file for Gromacs dynamics for $pdb"
   echo ""
   echo "ndx-maker.sh 0 OK" > ../arm.err
else
   echo " File $pdb not found!"
   echo ""
   echo "ndx-maker.sh 1 PDBNotFound" > ../arm.err
fi
#
# $moving=5 is to select the movable part for post-MD analysis (called by Analysis_MD.sh)
# Its value is changed by reading from Infos.dat the relevant data
#
if [[ $moving -eq 5 ]]; then
   realmove=$moving
   backbone=`grep "BackBoneMD" ../Infos.dat | awk '{ print $2 }'`
   modality=`grep "MDRelax" ../Infos.dat | awk '{ print $2 }'`
   case $modality in
   Hydrogen)
   moving=1
   ;;
   Cavity)
   moving=2
   ;;
   AllH)
   moving=3
   ;;
   AllSide)
   moving=4
   ;;
   esac
else
   realmove=0
fi

# Lines for atom to be fixed
# 
case $moving in
1)
if [[ $cavity == YES ]]; then
   selection="(all and not (($residline and hydrogen) or ((water within 4 of resname RET) and hydrogen)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2))"
else
   selection="all and not (hydrogen within $Hraggio of (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
fi
if [[ $risposta == y ]]; then
   selection="($selection) and not (resname RET and hydrogen)"
fi
;;
2) case $backbone in
   0) if [[ $cavity == YES ]]; then
         selection="all and not (($residline and sidechain) or (water within 4 of resname RET)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
      else
         selection="all and not ((same residue as (all within $radius of (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2))) and (sidechain or water)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
      fi
   ;;
   1) if [[ $cavity == YES ]]; then
         selection="all and not ($residline or (water within 4 of resname RET)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
      else
         selection="all and not ((same residue as (all within $radius of (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)))) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)" 
      fi
   ;;
   esac
   if [[ $risposta == y ]]; then
      selection="($selection) and not (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
   fi
;;
3)
selection="((all and not hydrogen) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) and not ions"
if [[ $risposta == y ]]; then
   selection="($selection) and not (resname RET and hydrogen)"
fi
;;
4)
selection="((all and not sidechain) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) and not (water or ions)"
if [[ $risposta == y ]]; then
   selection="($selection) and not (resname RET and sidechain)"
fi
;;
esac 
#
# Header for the ndx file
# The header is changed when $moving=5, i.e. when Analysis_MD.sh calls this script
# Also, when $moving=5 the selection has to be inverted
#
if [[ $realmove -eq 0 ]]; then
   echo "[ GroupDyna ]" > testa
else
   echo "[ Cavity ]" > testa
   truesel="all and not ($selection)"
   selection=$truesel
fi

# Storing the final selection for future use
#
echo $selection > cursel.log
 
# TCL script for VMD: open file, apply selection, save the serial numbers into a file
#
echo -e "mol new $pdb type gro" > ndxsel.tcl
echo -e "mol delrep 0 top" >> ndxsel.tcl
riga1="set babbeo [ atomselect top \"$selection\" ]"
echo -e "$riga1" >> ndxsel.tcl
echo -e 'set noah [$babbeo get serial]' >> ndxsel.tcl
riga3="set filename grodyn.dat"
echo -e "$riga3" >> ndxsel.tcl
echo -e 'set fileId [open $filename "w"]' >> ndxsel.tcl
echo -e 'puts -nonewline $fileId $noah' >> ndxsel.tcl
echo -e 'close $fileId' >> ndxsel.tcl
echo -e "exit" >> ndxsel.tcl
vmd -e ndxsel.tcl -dispdev text

# VIM script for putting grodyn.dat in the Gromacs ndx file format
#
echo ":set tw=75" > shiftline.vim
echo ":normal gqw" >> shiftline.vim
echo ":x" >> shiftline.vim
vim -es grodyn.dat < shiftline.vim
cat testa grodyn.dat >> $1.ndx
rm ndxsel.tcl shiftline.vim testa grodyn.dat

# If the cavity file is not provided and cavity is moved, the file is generated
#
if [[ $cavity != YES && $moving -eq 2 && $realmove -ne 5 ]]; then
   selection="all and not ($selection)"
cat > ndxsel.tcl << VMD1
mol new $pdb
mol delrep 0 top
set babbeo [ atomselect top "$selection" ]
set noah [\$babbeo get resid]
set nogo [\$babbeo get serial]
set filename residues.dat
set fileser serials.dat
set fileId [open \$filename "w"]
set file2Id [open \$fileser "w"]
puts -nonewline \$fileId \$noah
close \$fileId
puts -nonewline \$file2Id \$nogo
close \$file2Id
exit
VMD1
  vmd -e ndxsel.tcl -dispdev text
  sed -i 's/ /\n/g' residues.dat
  cat residues.dat | uniq > ../cavity
  rm residues.dat
  if [[ -z $cavity ]]; then
     ../update_infos.sh 0 "CavityFile" "YES" ../Infos.dat
  else
     ../update_infos.sh 1 "CavityFile" "YES" ../Infos.dat
  fi
fi

