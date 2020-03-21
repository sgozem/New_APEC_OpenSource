#!/bin/bash
#
checkrest=$1
chiave=$2
parametro=$3
infospath=$4
if [[ $checkrest -eq 0  ]]; then
   echo "$chiave $parametro" >> $infospath
else
   awk '{ if ( $1 == "'"$chiave"'" ) sub($2,"'"$parametro"'"); print $0 }' $infospath >> $infospath.vero
   mv $infospath.vero $infospath
fi

