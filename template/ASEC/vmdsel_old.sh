#!/bin/bash
#
#
# VMD selection: sidechains within 4 angstrom of QM region
#
pdbfile=$1.pdb
qmserials=`cat qmserials`
selection="((same residue as (all within 4 of serial $qmserials)) and (sidechain or water)) and not (serial $qmserials)"

# TCL script for VMD: open file, apply selection, save the serial numbers into a file
#
echo -e "mol new $pdbfile" > keysel.tcl
echo -e "mol delrep 0 top" >> keysel.tcl
riga1="set babbeo [ atomselect top \"$selection\" ]"
echo -e "$riga1" >> keysel.tcl
echo -e 'set noah [$babbeo get serial]' >> keysel.tcl
riga3="set filename qmmm.dat"
echo -e "$riga3" >> keysel.tcl
echo -e 'set fileId [open $filename "w"]' >> keysel.tcl
echo -e 'puts -nonewline $fileId $noah' >> keysel.tcl
echo -e 'close $fileId' >> keysel.tcl
echo -e "exit" >> keysel.tcl
vmd -e keysel.tcl -dispdev text
rm keysel.tcl
