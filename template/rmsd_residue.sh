#!/bin/bash
#
#
# This script generates the list of unique residues from a given selection 
#
pdb=$1
selection=$2
backb=$3
primo=$4
timeprod=$5
cat > ndxsel.tcl << VMD1
mol new $pdb.gro
mol delrep 0 top
set babbeo [ atomselect top "$selection" ]
set noah [\$babbeo get resid]
set filename residues.dat
set fileId [open \$filename "w"]
puts -nonewline \$fileId \$noah
close \$fileId
exit
VMD1
vmd -e ndxsel.tcl -dispdev text
sed -i 's/ /\n/g' residues.dat
cat residues.dat | uniq > listares
lista=( $( cat listares ) )
cat listares
sleep 1
rm residues.dat ndxsel.tcl listares
k=1
echo "set term png" > resplotter.gp
echo 'set output "'"group$k.png"'"' >> resplotter.gp
echo 'unset key' >> resplotter.gp
echo "set multiplot layout 2, 2" >> resplotter.gp
num=${#lista[@]}
num=$(($num-1))
for i in `seq 0 $num`; do
      residuo=${lista[$i]}
      if [[ $backb -eq 1 ]]; then
         selection="resid $residuo"
      else
         selection="resid $residuo and (sidechain or water)"
      fi
      cat > res-rmsd.tcl << MORO
mol new $pdb.gro
mol addfile $pdb.trr waitfor all
set steps [ molinfo 0 get numframes ]
set ref [ atomselect top "$selection" frame $primo ]
set compare [ atomselect top "$selection" ]
set file [ open "$residuo-rmsd.dat" "w" ]
for {set frame $primo} {\$frame < \$steps} {incr frame} {
    \$compare frame \$frame
    set rmsd [ measure rmsd \$compare \$ref ]
    puts \$file \$rmsd
}
close \$file
set beta [ atomselect top "resid $residuo and name CA OW" ] 
set nome [ \$beta get resname ]
set file2 [ open "nome$residuo.dat" "w" ]
puts \$file2 \$nome
close \$file2
exit      
MORO
    vmd -e res-rmsd.tcl -dispdev text
    resname=`cat nome$residuo.dat`
    rm nome$residuo.dat
    resto=$(($i % 4))
    if [[ $resto -eq 0 && $i -ne 0 ]]; then
       k=$(($k+1))
       echo 'unset multiplot' >> resplotter.gp
       echo 'set output "'"group$k.png"'"' >> resplotter.gp
       echo 'set multiplot layout 2, 2' >> resplotter.gp
    fi 
    echo 'set title "'"$resname$residuo"'"' >> resplotter.gp
    echo 'plot [0:'"$timeprod"'] [0:2] "'"$residuo-rmsd.dat"'" with lines' >> resplotter.gp
done
echo 'unset multiplot' >> resplotter.gp
gnuplot resplotter.gp
rm resplotter.gp res-rmsd.tcl
mkdir RMSD_Data Images
mv *-rmsd.dat RMSD_Data
mv group*.png Images 
