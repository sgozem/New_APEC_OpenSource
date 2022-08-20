#!/bin/bash
#
#
# VMD selection: sidechains within 4 angstrom of QM region or atoms moved during MD according to cavity residues file
#
pdbfile=$1.pdb
qmserials=`cat qmserials`
selection=$2
backb=$3
radius=$4
if [[ $selection == NO ]]; then
   if [[ $backb -eq 0 ]]; then
      selection="((same residue as (all within $radius of serial $qmserials)) and (sidechain or water)) and not (serial $qmserials)"
#      echo "VMDsel: $selection"
   else
      selection="((same residue as (all within $radius of serial $qmserials))) and not (serial $qmserials)"
#      echo "VMDsel: $selection"
   fi
else
#   if [[ $backb -eq 0 ]]; then
#      selection="($selection and (sidechain or water)) and not (serial $qmserials)"
#      echo "VMDsel: $selection"
#   else
#      selection="$selection and not (serial $qmserials)"
#      echo "VMDsel: $selection"
#   fi
   residline=`tr '\n' ' ' < ../cavity`
   if [[ $backb -eq 0 ]]; then
      backsel="and sidechain"
   fi
   acqua="or ((same residue as (all within 4 of (resname RET and not name N H CA HA C O CB CG CD HB1 HB2 HG1 HG2 HD1 HD2))) and water)"
   selection="(resid $residline $backsel) $acqua"
#   echo "VMDsel: $selection"
fi

# TCL script for VMD: open file, apply selection, save the serial numbers into a file
#
#echo -e "mol new $pdbfile" > keysel.tcl
#echo -e "mol delrep 0 top" >> keysel.tcl
#riga1="set babbeo [ atomselect top \"$selection\" ]"
#echo -e "$riga1" >> keysel.tcl
#echo -e 'set noah [$babbeo get serial]' >> keysel.tcl
#riga3="set filename qmmm.dat"
#echo -e "$riga3" >> keysel.tcl
#echo -e 'set fileId [open $filename "w"]' >> keysel.tcl
#echo -e 'puts -nonewline $fileId $noah' >> keysel.tcl
#echo -e 'close $fileId' >> keysel.tcl
#echo -e "exit" >> keysel.tcl
cat > keysel.tcl << MOROK
mol new $pdbfile
mol delrep 0 top
set babbeo [ atomselect top "$selection" ]
set noah [ \$babbeo get serial ]
set filename qmmm.dat
set fileId [open \$filename "w"]
puts -nonewline \$fileId \$noah
close \$fileId
exit
MOROK
vmd -e keysel.tcl -dispdev text
rm keysel.tcl
