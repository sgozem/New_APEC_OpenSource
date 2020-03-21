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

cavity=`grep "CavityFile" ../Infos.dat | awk '{ print $2 }'`
radius=`grep "RadiusMD" ../Infos.dat | awk '{ print $2 }'`

Hraggio=7

#fi
#if [[ -z $radius ]]; then
#   radius=4
#fi
#if [[ -z $Hraggio ]]; then
#   Hraggio=7
#fi

if [[ $cavity == YES ]]; then
   waters=`grep "Waters" ../Infos.dat | awk '{ print $2 }'`
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
   echo " Sidechains within "cavity" will be relaxed"
   echo ""
else
   echo " Sidechains within $radius will be relaxed"
   echo ""  
fi

#
# Checking for PDB file existence
#
if [ -f $pdb ]; then
   echo " Creating the ndx file for Gromacs dynamics for $pdb"
   echo ""
   echo "ndx-maker OK" > log
else
   echo " File $pdb not found!"
   echo ""
   echo "ndx-maker failed" > log
fi

# Lines for atom to be fixed
# 
case $moving in
1)
selection="all and not (hydrogen within $Hraggio of (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)"
;;
2) case $backbone in
   0) if [[ $cavity == YES ]]; then
         selection="all and not (($residline and (sidechain or (resname ACI ACH))) or (water within $waters of (resname RET and not name N H CA HA C O))) or (resname RET and not name N H CA HA C O)"
      else
         selection="all and not ((same residue as (all within $radius of (resname RET and not name N H CA HA C O))) and (sidechain or water or (resname ACI ACH))) or (resname RET and not name N H CA HA C O)"
      fi
   ;;
# We need maintain fixed the C and N of the previous and next backbone aminoacid in order to be able of
# optimizing (MM) the C(B) of the RET because of the dihidral angles. That is the reason of the:
# ((all within 2.0 of (resname RET and name C)) and name N) or ((all within 2.0 of (resname RET and name N)) and name C)

   1) if [[ $cavity == YES ]]; then
         selection="all and not (($residline) or (water within $waters of (resname RET and not name N H CA HA C O))) or (resname RET) or (((all within 2.0 of (resname RET and name C)) and name N) or ((all within 2.0 of (resname RET and name N)) and name C))"
      else
         selection="all and not (same residue as (all within $radius of (resname RET and not name N H CA HA C O))) or (resname RET) or ((all within 2.0 of (resname RET and name C)) and name N) or ((all within 2.0 of (resname RET and name N)) and name C)"
      fi 
   ;;
   esac
;;
3)
selection="((all and not hydrogen) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) and not ions"
;;
4)
selection="((all and not sidechain) or (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2)) and not (water or ions)"
;;
5) case $backbone in
   0) selection="backbone or resname CHR" 
   ;;
   1) selection="resname CHR"
   ;;
   esac
;;
esac 
 
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

cp grodyn.dat list

num=`awk '{print NF}' list`
echo " $num" > list2
echo " " >> list2
for i in $(eval echo "{1..$num}")
do
  awk -v j=$i '{ print $j }' list >> list2
done

mv list2 list_gro
rm list

# Header for the ndx file
#
echo "[ GroupDyna ]" > testa

# VIM script for putting grodyn.dat in the Gromacs ndx file format
#
echo ":set tw=75" > shiftline.vim
echo ":normal gqw" >> shiftline.vim
echo ":x" >> shiftline.vim
vim -es grodyn.dat < shiftline.vim
cat testa grodyn.dat >> $1.ndx
rm ndxsel.tcl shiftline.vim testa grodyn.dat

#selec=`head -n1 list | awk '{ print $0 }'`

if [[ $cavity != YES ]]; then
selection="all and not ($selection)"
#selection="all and not (serial $selec)"
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
#  if [[ -z $cavity ]]; then
     ../update_infos.sh "CavityFile" "YES" ../Infos.dat
#  else
#     ../update_infos.sh 1 "CavityFile" "YES" ../Infos.dat
#  fi
fi
