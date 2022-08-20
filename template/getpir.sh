#!/bin/bash
#
pdb="PROGETTO.pdb"
quanti=`grep 'ATOM ' $pdb | tail -n 1 | awk '{ print $6 }'`
flag=0
# YOE bug initiate resall
for i in $(seq 1 $quanti); do
    residuo=`awk '{ if ( $6 == "'"$i"'" ) print $4 }' $pdb | uniq`
    if [[ ! -z $residuo ]]; then
       ./aminoconv.sh $residuo $residuo 1
       reslett1=`cat wtres`
       resall=$resall$reslett1
       flag=0
    else
       if [[ $flag -eq 0 ]]; then
          resall=$resall"/"
          flag=1
       fi
    fi
done
rm wtres sidemut
echo $resall > PROGETTO.pir

