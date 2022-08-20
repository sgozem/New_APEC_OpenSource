#!/bin/bash

echo ""
echo ""
echo " What is the name of the pdb file to convert (omit the extension)"
echo ""
read file

cp $file.pdb new_$file.pdb

sed -i "/CONECT/d" new_$file.pdb
sed -i "/ H /d;/ 1H. /d; / 2H. /d; / 3H. /d; / HA/d; / HB/d; / HG/d; / HD/d; / HE/d; / HH/d; / HZ/d; / HT/d" new_$file.pdb
sed -i "/ OXT /d" new_$file.pdb

# Retrieving the retinal-bound lysine number (From Federico, New.sh)
#
x=`grep 'C15 RET' new_$file.pdb | awk '{ print $7 }'`
y=`grep 'C15 RET' new_$file.pdb | awk '{ print $8 }'`
z=`grep 'C15 RET' new_$file.pdb | awk '{ print $9 }'`
lysnum=`grep "LYS" new_$file.pdb | grep " NZ " | awk '{ xdiff = ($7 - "'"$x"'"); ydiff = ($8 - "'"$y"'"); zdiff = ($9 - "'"$z"'"); if ( xdiff < 0 ) xdiff=-xdiff; if ( ydiff < 0 ) ydiff=-ydiff; if ( zdiff < 0 ) zdiff=-zdiff;  if ( xdiff < 1.5 && ydiff < 1.5 && zdiff < 1.5 ) print $6 }'`
echo " The number of the retinal-bound lysine is $lysnum"
echo " If it is ok type y, otherwise type the right number"
read risposta
echo ""
while [[ $risposta != [0-9][0-9][0-9] && $risposta != y ]]; do
   echo " Wrong characters! Please type a 3-digit number"
   read risposta
   echo ""
done
if [[ $risposta == [0-9][0-9][0-9] ]]; then
   lysnum=$risposta
fi

chain=`grep " NZ  LYS . $lysnum" new_$file.pdb | awk '{ print $5 }'`

sed -i "s/ RET Z   1 / RET Z 999 /g" new_$file.pdb
sed -i "/RET/s/ATOM  /HETATM/; /RET/s/ Z / $chain /" new_$file.pdb
sed -i "s/ OT  WAT  / O   HOH $chain/" new_$file.pdb

answer=b
while [[ $answer != "y" && $answer != "n" ]]; do
   echo ""
   echo " Do you want to keep the CL or NA ions in the pdb file?"
   echo ""
   echo ""
   read answer
done
if [[ $answer == "n" ]]; then
   sed -i "/NA    NA /d" new_$file.pdb
   sed -i "/CL    CL /d" new_$file.pdb
fi

for i in {21..48}; do
   sed -i "/ H$i /d" new_$file.pdb
done

sed -i "/ C1  RET A/iTER" new_$file.pdb

line=`grep -n "C20 RET $chain" new_$file.pdb | awk '{ print $1 }' FS=":"`

head -n $(($line+20)) new_$file.pdb > temp
mv temp new_$file.pdb


