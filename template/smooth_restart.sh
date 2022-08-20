#!/bin/bash
#
folder=$1
message=$2
step=$3
echo " Folder $folder was found! $message"
read conferma
echo ""
if [[ $conferma == y ]]; then
   i=1
   while [[ -d $i.$folder ]]; do
         i=$(($i+1))
   done
   mv $folder $i.$folder
   cp Infos.dat $i.Infos.dat
   echo "ReStRun_$step $i" >> Infos.dat
   cp arm.err $i.arm.err
else
   echo " No further operation! Aborting..."
   echo ""
   mv Infos.dat no.Infos.dat
fi

