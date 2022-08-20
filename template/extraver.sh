#!/bin/bash
#
Project=$1
primo=$2
residfile=$3
backb=$4
residline=`tr '\n' ' ' < $residfile`
if [[ $backb -eq 0 ]]; then
   selection="resid $residline and (sidechain or water)" 
else
   selection="resid $residline"
fi
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`
echo 0 > choices.txt
$gropath/gmx g_rmsf -f $Project.trr -s $Project.tpr -n $Project.ndx -ox < choices.txt
rm choices.txt
cat > extr-aver.tcl << KUMA
mol new xaver.pdb
mol new $Project.gro
mol addfile $Project.trr waitfor all
set ref [ atomselect 0 "$selection" ]
set compare [ atomselect top "$selection" ]
set steps [ molinfo top get numframes ]
set file [ open "lowrmsd.dat" "w" ]
set min_rmsd 100
for {set frame $primo} {\$frame < \$steps} {incr frame} {
    \$compare frame \$frame
    set rmsd [ measure rmsd \$compare \$ref ]
    puts \$file \$rmsd
    if { \$rmsd <= \$min_rmsd } {
       set min_rmsd \$rmsd
       set min_frame \$frame
    }
}
set min_frame [ expr \$min_frame - 1 ]
puts \$file \$min_frame
close \$file
exit
KUMA
vmd -e extr-aver.tcl -dispdev text
rm extr-aver.tcl
