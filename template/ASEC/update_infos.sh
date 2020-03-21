#!/bin/bash
#
parameter=$1
value=$2
file=$3

key=`grep "$parameter" $file | awk '{ print $1 }'`

if [[ $key == $parameter ]]; then
   awk '{ if ( $1 == "'"$parameter"'" ) sub($2,"'"$value"'"); print $0 }' $file >> Infos_ok.dat
   mv Infos_ok.dat $file
else
   echo "$parameter $value" >> $file
fi

