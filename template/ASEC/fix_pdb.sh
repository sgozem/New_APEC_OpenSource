#!/bin/bash
#

echo ""
echo " What is the name of the pdb file obtained after running the Standard ARM, (omit the extension .pdb)"
read Project
echo ""

#echo ""
#echo " What is the full path of the ARM protocol (i.e. .../template)"
#read camino
#echo ""

#Define the ARM template folder here:
camino=/users/PCS0202/bgs0361/bin/ARM_protocol/template

######
# Retrieving the retinal-bound lysine number
########## From Federico at New.sh #########
######

x=`grep 'C15 RET' ${Project}.pdb | awk '{ print $7 }'`
y=`grep 'C15 RET' ${Project}.pdb | awk '{ print $8 }'`
z=`grep 'C15 RET' ${Project}.pdb | awk '{ print $9 }'`
Lys=`grep "LYS" ${Project}.pdb | grep " NZ " | awk '{ xdiff = ($7 - "'"$x"'"); ydiff = ($8 - "'"$y"'"); zdiff = ($9 - "'"$z"'"); if ( xdiff < 0 ) xdiff=-xdiff; if ( ydiff < 0 ) ydiff=-ydiff; if ( zdiff < 0 ) zdiff=-zdiff;  if ( xdiff < 1.5 && ydiff < 1.5 && zdiff < 1.5 ) print $6 }'`
echo " The number of the retinal-bound lysine is $Lys"
echo " If it is ok type y, otherwise type the right number"
read risposta
echo ""
while [[ $risposta != [0-9][0-9][0-9] && $risposta != y ]]; do
   echo " Wrong characters! Please type a 3-digit number"
   read risposta
   echo ""
done
if [[ $risposta == [0-9][0-9][0-9] ]]; then
   Lys=$risposta
fi

################################

if [ $Lys -lt 10 ]; then
   Chain=`grep " HZ1 LYS.....$Lys " $Project.pdb | awk '{ print $5 }'`
fi
if [ $Lys -ge 10 ] && [ $Lys -lt 100 ]; then
   Chain=`grep " HZ1 LYS....$Lys " $Project.pdb | awk '{ print $5 }'`
fi
if [ $Lys -ge 100 ] && [ $Lys -lt 1000 ]; then
   Chain=`grep " HZ1 LYS...$Lys " $Project.pdb | awk '{ print $5 }'`
fi

grep " RET Z   1 " $Project.pdb > RETI
#sed -i "s/ RET Z   1 / RET $Chain $Lys /g" RETI  
sed "/ RET Z   1 /d" $Project.pdb > ${Project}_fixed.pdb

grep " ACI " ${Project}_fixed.pdb > CLORI
sed -i "/ ACI /d" ${Project}_fixed.pdb
grep " ACH " ${Project}_fixed.pdb >> CLORI
sed -i "/ ACH /d" ${Project}_fixed.pdb

grep " CL " ${Project}_fixed.pdb >> CLORI
grep " NA " ${Project}_fixed.pdb >> CLORI
sed -i "/ CL /d" ${Project}_fixed.pdb
sed -i "/ NA /d" ${Project}_fixed.pdb
sed -i "/CONECT/d" ${Project}_fixed.pdb
sed -i "/END/d" ${Project}_fixed.pdb
sed -i "/HEADER/d" ${Project}_fixed.pdb
sed -i "/COMPND/d" ${Project}_fixed.pdb
sed -i "/SOURCE/d" ${Project}_fixed.pdb
grep " WAT " ${Project}_fixed.pdb > WATER
sed -i "/ WAT /d" ${Project}_fixed.pdb

if [ $Lys -lt 10 ]; then
   sed -i "s/ RET Z   1 / RET $Chain   $Lys /g" RETI
   sed -i "s/LYS "$Chain"   "$Lys"/RET "$Chain" "$Lys"/g" ${Project}_fixed.pdb
fi
if [ $Lys -ge 10 ] && [ $Lys -lt 100 ]; then
   sed -i "s/ RET Z   1 / RET $Chain  $Lys /g" RETI
   sed -i "s/LYS "$Chain"  "$Lys"/RET "$Chain" "$Lys"/g" ${Project}_fixed.pdb
fi
if [ $Lys -ge 100 ] && [ $Lys -lt 1000 ]; then
   sed -i "s/ RET Z   1 / RET $Chain $Lys /g" RETI
   sed -i "s/LYS "$Chain" "$Lys"/RET "$Chain" "$Lys"/g" ${Project}_fixed.pdb
fi

#sed -i "s/LYS "$Chain" "$Lys"/RET "$Chain" "$Lys"/g" ${Project}_fixed.pdb

LastRet=`grep -n " HZ1 RET " ${Project}_fixed.pdb | awk '{ print $1 }' FS=":"`
head -n $LastRet ${Project}_fixed.pdb > temp1
cat temp1 RETI > temp2
sed -e "1,$LastRet d" ${Project}_fixed.pdb > temp3

wat=`grep -c " OT  WAT " WATER`
if [[ -f WATER2 ]]; then
   rm WATER2
fi

for i in $(eval echo "{1..$wat}")
do
  head -n $((3*$i-2)) WATER | tail -n1 > temp4
  sed -i "s/ OT  WAT / OW  HOH /g" temp4
  head -n1 temp4 >> WATER2
  
  head -n $(((3*$i-2)+1)) WATER | tail -n1 > temp4
  sed -i "s/ HT  WAT / HW1 HOH /g" temp4
  head -n1 temp4 >> WATER2
 
  head -n $(((3*$i-2)+2)) WATER | tail -n1 > temp4
  sed -i "s/ HT  WAT / HW2 HOH /g" temp4
  head -n1 temp4 >> WATER2
done

cat temp2 temp3 > temp5
echo "TER" >> temp5
cat temp5 WATER2 CLORI > ${Project}_fixed.pdb
rm temp* RETI CLORI WATER WATER2

#####################

sed -i "s/ HIP / HIS /g" ${Project}_fixed.pdb
sed -i "s/ HID / HIS /g" ${Project}_fixed.pdb
sed -i "s/ HIE / HIS /g" ${Project}_fixed.pdb

sed -i "s/ HB2 MET / HB1 MET /g" ${Project}_fixed.pdb
sed -i "s/ HB2 ASN / HB1 ASN /g" ${Project}_fixed.pdb
sed -i "s/ HB2 GLU / HB1 GLU /g" ${Project}_fixed.pdb
sed -i "s/ HB2 PRO / HB1 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HB2 PHE / HB1 PHE /g" ${Project}_fixed.pdb
sed -i "s/ HB2 TYR / HB1 TYR /g" ${Project}_fixed.pdb
sed -i "s/ HB2 SER / HB1 SER /g" ${Project}_fixed.pdb
sed -i "s/ HB2 LYS / HB1 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HB2 ARG / HB1 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HB2 GLN / HB1 GLN /g" ${Project}_fixed.pdb
sed -i "s/ HB2 LEU / HB1 LEU /g" ${Project}_fixed.pdb
sed -i "s/ HB2 TRP / HB1 TRP /g" ${Project}_fixed.pdb
sed -i "s/ HB2 HIS / HB1 HIS /g" ${Project}_fixed.pdb
sed -i "s/ HB2 CYS / HB1 CYS /g" ${Project}_fixed.pdb
sed -i "s/ HB2 ASP / HB1 ASP /g" ${Project}_fixed.pdb
sed -i "s/ HB2 RET / HB1 RET /g" ${Project}_fixed.pdb


sed -i "s/ HB3 MET / HB2 MET /g" ${Project}_fixed.pdb
sed -i "s/ HB3 ASN / HB2 ASN /g" ${Project}_fixed.pdb
sed -i "s/ HB3 GLU / HB2 GLU /g" ${Project}_fixed.pdb
sed -i "s/ HB3 PRO / HB2 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HB3 PHE / HB2 PHE /g" ${Project}_fixed.pdb
sed -i "s/ HB3 TYR / HB2 TYR /g" ${Project}_fixed.pdb
sed -i "s/ HB3 SER / HB2 SER /g" ${Project}_fixed.pdb
sed -i "s/ HB3 LYS / HB2 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HB3 ARG / HB2 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HB3 GLN / HB2 GLN /g" ${Project}_fixed.pdb
sed -i "s/ HB3 LEU / HB2 LEU /g" ${Project}_fixed.pdb
sed -i "s/ HB3 TRP / HB2 TRP /g" ${Project}_fixed.pdb
sed -i "s/ HB3 HIS / HB2 HIS /g" ${Project}_fixed.pdb
sed -i "s/ HB3 CYS / HB2 CYS /g" ${Project}_fixed.pdb
sed -i "s/ HB3 ASP / HB2 ASP /g" ${Project}_fixed.pdb
sed -i "s/ HB3 RET / HB2 RET /g" ${Project}_fixed.pdb

###########

sed -i "s/ HG2 MET / HG1 MET /g" ${Project}_fixed.pdb
sed -i "s/ HG2 PRO / HG1 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HG2 GLU / HG1 GLU /g" ${Project}_fixed.pdb
sed -i "s/ HG2 ARG / HG1 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HG2 LYS / HG1 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HG2 GLN / HG1 GLN /g" ${Project}_fixed.pdb
sed -i "s/ HG2 RET / HG1 RET /g" ${Project}_fixed.pdb

sed -i "s/ HG3 MET / HG2 MET /g" ${Project}_fixed.pdb
sed -i "s/ HG3 PRO / HG2 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HG3 GLU / HG2 GLU /g" ${Project}_fixed.pdb
sed -i "s/ HG3 ARG / HG2 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HG3 LYS / HG2 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HG3 GLN / HG2 GLN /g" ${Project}_fixed.pdb
sed -i "s/ HG3 RET / HG2 RET /g" ${Project}_fixed.pdb

###########

sed -i "s/ HD2 PRO / HD1 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HD2 LYS / HD1 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HD2 ARG / HD1 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HD2 RET / HD1 RET /g" ${Project}_fixed.pdb

sed -i "s/ HD3 PRO / HD2 PRO /g" ${Project}_fixed.pdb
sed -i "s/ HD3 LYS / HD2 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HD3 ARG / HD2 ARG /g" ${Project}_fixed.pdb
sed -i "s/ HD3 RET / HD2 RET /g" ${Project}_fixed.pdb

###########

sed -i "s/ HE2 LYS / HE1 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HE2 RET / HE1 RET /g" ${Project}_fixed.pdb

sed -i "s/ HE3 LYS / HE2 LYS /g" ${Project}_fixed.pdb
sed -i "s/ HE3 RET / HE2 RET /g" ${Project}_fixed.pdb


sed -i "s/ 1HA  GLY /  HA1 GLY /g" ${Project}_fixed.pdb
sed -i "s/ HA3 GLY / HA2 GLY /g" ${Project}_fixed.pdb


sed -i "s/ HG12 ILE / HG11 ILE /g" ${Project}_fixed.pdb
sed -i "s/ HG13 ILE / HG12 ILE /g" ${Project}_fixed.pdb

###########
###########
#RET
###########
###########

sed -i "s/ C10 RET / CF  RET /" ${Project}_fixed.pdb
sed -i "s/ C11 RET / CI  RET /" ${Project}_fixed.pdb
sed -i "s/ C12 RET / CJ  RET /" ${Project}_fixed.pdb
sed -i "s/ C14 RET / CK  RET /" ${Project}_fixed.pdb
sed -i "s/ C15 RET / CM  RET /" ${Project}_fixed.pdb
sed -i "s/ C16 RET / CO  RET /" ${Project}_fixed.pdb
sed -i "s/ C17 RET / CP  RET /" ${Project}_fixed.pdb
sed -i "s/ C18 RET / CQ  RET /" ${Project}_fixed.pdb
sed -i "s/ C19 RET / CU  RET /" ${Project}_fixed.pdb
sed -i "s/ C20 RET / CX  RET /" ${Project}_fixed.pdb

sed -i "s/ H21 RET /HC21 RET /" ${Project}_fixed.pdb
sed -i "s/ H22 RET /HC22 RET /" ${Project}_fixed.pdb
sed -i "s/ H23 RET /HC31 RET /" ${Project}_fixed.pdb
sed -i "s/ H24 RET /HC32 RET /" ${Project}_fixed.pdb
sed -i "s/ H25 RET /HC41 RET /" ${Project}_fixed.pdb
sed -i "s/ H26 RET /HC42 RET /" ${Project}_fixed.pdb
sed -i "s/ H27 RET /HC71 RET /" ${Project}_fixed.pdb
sed -i "s/ H28 RET /HC81 RET /" ${Project}_fixed.pdb
sed -i "s/ H29 RET /HCF1 RET /" ${Project}_fixed.pdb
sed -i "s/ H30 RET /HCI1 RET /" ${Project}_fixed.pdb
sed -i "s/ H31 RET /HCJ1 RET /" ${Project}_fixed.pdb
sed -i "s/ H32 RET /HCK1 RET /" ${Project}_fixed.pdb
sed -i "s/ H33 RET /HCM1 RET /" ${Project}_fixed.pdb
sed -i "s/ H34 RET /HCO1 RET /" ${Project}_fixed.pdb
sed -i "s/ H35 RET /HCO2 RET /" ${Project}_fixed.pdb
sed -i "s/ H36 RET /HCO3 RET /" ${Project}_fixed.pdb
sed -i "s/ H37 RET /HCP1 RET /" ${Project}_fixed.pdb
sed -i "s/ H38 RET /HCP2 RET /" ${Project}_fixed.pdb
sed -i "s/ H39 RET /HCP3 RET /" ${Project}_fixed.pdb
sed -i "s/ H40 RET /HCQ1 RET /" ${Project}_fixed.pdb
sed -i "s/ H41 RET /HCQ2 RET /" ${Project}_fixed.pdb
sed -i "s/ H42 RET /HCQ3 RET /" ${Project}_fixed.pdb
sed -i "s/ H43 RET /HCU1 RET /" ${Project}_fixed.pdb
sed -i "s/ H44 RET /HCU2 RET /" ${Project}_fixed.pdb
sed -i "s/ H45 RET /HCU3 RET /" ${Project}_fixed.pdb
sed -i "s/ H46 RET /HCX1 RET /" ${Project}_fixed.pdb
sed -i "s/ H47 RET /HCX2 RET /" ${Project}_fixed.pdb
sed -i "s/ H48 RET /HCX3 RET /" ${Project}_fixed.pdb

grep " HIS " ${Project}_fixed.pdb > temp1

ii=`wc -l temp1 | awk '{ print $1 }'`
if [[ $ii -gt 0 ]]; then
   cont=1

   n=`head -n1 temp1 | awk '{ print $6 }'`
   hischain=`head -n1 temp1 | awk '{ print $5 }'`
   echo "$n $hischain"> temp2

   for i in $(eval echo "{2..$ii}")
   do
    k=`head -n $i temp1 | tail -n1 | awk '{ print $6 }'`
    hischain=`head -n $i temp1 | tail -n1 | awk '{ print $5 }'`
    if [[ $k != $n ]]; then
     echo "$k $hischain" >> temp2
     n=$k
    fi
   done
   rm temp1

   ii=`wc -l temp2 | awk '{ print $1 }'`
   i=1
   for i in $(eval echo "{1..$ii}")
   do
     his=`head -n $i temp2 | tail -n1 | awk '{ print $1 }'`
     hischain=`head -n $i temp2 | tail -n1 | awk '{ print $2 }'`
     echo ""
     echo " Please define the protonation state of the following $ii HISs:"
     echo " HIS $his"
     echo ""
     if [ $his -lt 10 ]; then
       grep " HIS.....$his " ${Project}_fixed.pdb
     fi
     if [ $his -ge 10 ] && [ $his -lt 100 ]; then
       grep " HIS....$his " ${Project}_fixed.pdb
     fi
     if [ $his -ge 100 ] && [ $his -lt 1000 ]; then
       grep " HIS...$his " ${Project}_fixed.pdb
     fi
 
     answer=0
     while  [[ $answer -ne 1 && $answer -ne 2 && $answer -ne 3 ]]; do
           echo " Select one of the following protonation states if the case"
           echo ""
           echo " 1) HIP (18 atoms)"
           echo " 2) HID (17 atoms, with HD1 and HD2)"
           echo " 3) HIE (17 atoms, with HE1 and HE2)"
           echo ""
           read answer
     done
 
     if [[ $answer == 1 ]]; then
        if [ $his -lt 10 ]; then
          sed -i "s/ HIS $hischain   $his / HIP $hischain   $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 10 ] && [ $his -lt 100 ]; then
          sed -i "s/ HIS $hischain  $his / HIP $hischain  $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 100 ] && [ $his -lt 1000 ]; then
          sed -i "s/ HIS $hischain $his / HIP $hischain $his /g" ${Project}_fixed.pdb
        fi
     fi
     if [[ $answer == 2 ]]; then
        if [ $his -lt 10 ]; then
          sed -i "s/ HIS $hischain   $his / HID $hischain   $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 10 ] && [ $his -lt 100 ]; then
          sed -i "s/ HIS $hischain  $his / HID $hischain  $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 100 ] && [ $his -lt 1000 ]; then
          sed -i "s/ HIS $hischain $his / HID $hischain $his /g" ${Project}_fixed.pdb
        fi
     fi
     if [[ $answer == 3 ]]; then
        if [ $his -lt 10 ]; then
          sed -i "s/ HIS $hischain   $his / HIE $hischain   $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 10 ] && [ $his -lt 100 ]; then
          sed -i "s/ HIS $hischain  $his / HIE $hischain  $his /g" ${Project}_fixed.pdb
        fi
        if [ $his -ge 100 ] && [ $his -lt 1000 ]; then
          sed -i "s/ HIS $hischain $his / HIE $hischain $his /g" ${Project}_fixed.pdb
        fi
     fi
   done
   rm temp2
else
   rm temp1
   echo ""
   echo " There is no HIS sidechains in the PDB file"
   echo ""
fi
echo ""
echo " wait ..."
echo ""

#
# Inserting TER between two chains
#
atoms=`grep -c "ATOM " ${Project}_fixed.pdb`
chain="A"
cont=0
cp ${Project}_fixed.pdb temmp
for i in $(eval echo "{1..$atoms}")
do
   newchain=`head -n $i ${Project}_fixed.pdb | tail -n1 | awk '{ print $5 }'`
   if [[ $newchain != $chain ]]; then
      chain=$newchain 
      sed "$(($i+$cont))iTER" temmp > temmp2
      cont=$cont+1
      mv temmp2 temmp
   fi
done

mv temmp ${Project}_fixed.pdb

mv $Project.pdb ${Project}_QMMM.pdb
mv ${Project}_fixed.pdb $Project.pdb

cp $camino/ASEC/New_mod.sh .
sed -i "s|camino|$camino|g" New_mod.sh
sed -i "s|Projecto|$Project|" New_mod.sh
sed -i "s|lisina|$Lys|" New_mod.sh

echo ""
echo " Continues with the standard New_mod.sh"
echo "" 
###########
###########
# HIS
###########
###########
#sed -i "s/ HIS A 100 / HIP A 100 /g" ${Project}_fixed.pdb
#sed -i "s/ HIS A 152 / HIP A 152 /g" ${Project}_fixed.pdb
#sed -i "s/ HIS A 195 / HIP A 195 /g" ${Project}_fixed.pdb
#sed -i "s/ HIS A 211 / HID A 211 /g" ${Project}_fixed.pdb
#sed -i "s/ HIS A 278 / HIP A 278 /g" ${Project}_fixed.pdb



