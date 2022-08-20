#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
dowser=`grep "Dowser" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
chr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`
chargechr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
fmnfad=`grep "Tail" Infos.dat | awk '{ print $2 }'`

#########################################################################
# Please modify here any time new chromophore with new atom types is used
# NC, C, O, ... correspond to the gromacs atom types, while =1022, =3, ...
# correspond to the corresponding atom type in TINKER format.
#########################################################################
declare -A force
force=(["NC"]=1022 ["C"]=3 ["O"]=5 ["CA"]=115 ["CT"]=2 ["H"]=4 ["HC"]=14 ["Nstar"]=1017 ["H1"]=6 ["OH"]=63 ["HO"]=64 ["OS"]=1239 ["O2"]=1236 ["P"]=1235 ["CK"]=1021 ["NB"]=193 ["CB"]=149 ["N2"]=299 ["CQ"]=1023 ["H5"]=175)

#
# ESPF parametrization of the chromophore
#
mkdir ESPF_charges
cp Chromophore/${chromophore}.xyz ESPF_charges
cp Templates/$Project-tk.xyz ESPF_charges
cp Templates/$prm.prm ESPF_charges

cd ESPF_charges
#
# Defining the xyz and adding the connectivities to the Chromophore
#

numchr=`head -n1 ${chromophore}.xyz | awk '{ print $1 }'`
last=`head -n1 $Project-tk.xyz | awk '{ print $1 }'`

echo "$(($last+$numchr))" > top
head -n $(($last+1)) $Project-tk.xyz | tail -n $last >> top
cp ${chromophore}.xyz chromo.xyz
sed -i "s/\*/star/g" chromo.xyz
for i in $(eval echo "{1..$numchr}"); do
   att=`head -n $(($i+2)) chromo.xyz | tail -n1 | awk '{ print $5 }'`
   att=${force[$att]}
   x=`head -n $(($i+2)) ${chromophore}.xyz | tail -n1 | awk '{ print $2 }'`
#   x=$(echo "scale=2; ($x*10.0)" | bc)
   y=`head -n $(($i+2)) ${chromophore}.xyz | tail -n1 | awk '{ print $3 }'`
#   y=$(echo "scale=2; ($y*10.0)" | bc)
   z=`head -n $(($i+2)) ${chromophore}.xyz | tail -n1 | awk '{ print $4 }'`
#   z=$(echo "scale=2; ($z*10.0)" | bc)
   head -n $(($i+2)) ${chromophore}.xyz | tail -n1 | awk -v var1="$(($last+$i))" -v var2="$att" -v x=$x -v y=$y -v z=$z '{ print "  "var1"  "$1"   "x"       "y"       "z"     "var2 }' >> top
done
echo "" >> top

mv $Project-tk.xyz $Project-tk_noCHR.xyz
mv top $Project-tk.xyz

$tinkerdir/xyzedit $Project-tk.xyz << EOF
../$prm
7
EOF

rm chromo.xyz

echo "########################################################################"
echo ""
echo " It is normal the code to show 0 conectivities of the chromophore above."
echo " Indeed, the conectivities are being added now based on distace"
echo ""
echo "########################################################################"

mv $Project-tk.xyz_2 $Project-tk.xyz

cp $templatedir/ASEC/charges_$fmnfad template_charges

echo "parameters $prm" > $Project.key
echo "parameters $prm" > ${Project}_nolink.key
echo "" >> $Project.key
echo "" >> ${Project}_nolink.key
xx=0
fx=0
final=`head -n1 $Project-tk.xyz | awk '{ print $1 }'`
init=$(($final-$numchr))
rm -f mm
rm -f qm
rm -f chargemm
rm -f chargexx
rm -f chargefx
rm -f chargefxx
touch mm
touch qm
for i in $(eval echo "{1..$numchr}"); do
   num=`head -n $(($init+1+$i)) $Project-tk.xyz | tail -n1 | awk '{ print $1 }'`
   charge=`head -n $i template_charges | tail -n1 | awk '{ print $1 }'`
   atmtype=`head -n $(($i+2)) $chromophore.xyz | tail -n1 | awk '{ print $6 }'`
   if [[ $atmtype == "MM" ]]; then
      echo "MM $num" >> mm
      echo "CHARGE  -$num   $charge" >> chargemm
   fi
   if [[ $atmtype == "QM" ]]; then
      echo "QM $num" >> qm
   fi
   if [[ $atmtype == "LQ" ]]; then
      echo "QM $num" >> qm
   fi
   if [[ $atmtype == "LM" ]]; then
      echo "MM $num" >> mm
      echo "CHARGE  -$num   0.0000" >> chargemm
      LM=$num
   fi
#
#  1 for atoms of the tail as ASEC points
#  0 for the fixed atoms of the tail
#
   if [[ $atmtype == "XX" ]]; then
      echo "CHARGE  $num   $charge    1"   >> chargefxx
      xx=$(($xx+1))
   fi
   if [[ $atmtype == "FX" ]]; then
      echo "CHARGE  $num   $charge    0" >> chargefxx
      fx=$(($fx+1))
   fi
done
if [[ $xx -gt 0 ]]; then
   correct=a
   while [[ $correct != "y" && $correct != "n" ]]; do
      echo ""
      echo " *************************************************************************"
      echo ""
      echo " It seems that the chromophore has a long tail, which just the LM atom"
      echo " will be considered in the QMMM calculations as MM atoms. The rest of the"
      echo " tail will be considered as ASEC points."
      echo ""
      echo " Is it correct? (y/n)"
      echo ""
      echo " *************************************************************************"
      read correct
   done
   if [[ $correct == "n" ]]; then
      echo " "
      echo " Please redifine the QM, MM, QL, ML, and XX atom types in the initial"
      echo " chromophore xyz file."
      echo " Aborting ..."
      exit 0
   fi
fi

if [[ $xx -gt 0 ]]; then
   lines=`wc -l chargefxx | awk '{ print $1 }'`

   rm -f col
   for i in $(eval echo "{1..$lines}"); do
      tipo=`head -n $i chargefxx | tail -n1 | awk '{ print $4 }'`
      echo "  $tipo" >> col
   done
fi

if [[ $xx -gt 0 ]]; then
   paste chargefxx col > colu
   mv colu chargefxx
fi

#
# continuing generating the key file
#
echo "QMMM $(($numchr-$xx-$fx+1))" >> $Project.key
echo "QMMM $(($numchr-$xx-$fx))" >> ${Project}_nolink.key

cat mm >> $Project.key
cat mm >> ${Project}_nolink.key
cat qm >> $Project.key
cat qm >> ${Project}_nolink.key
echo "LA $(($final+1))" >> $Project.key

#echo "QMMM-ELECTROSTATICS ESPF" >> $Project.key
echo "QMMM-microiteration ON" >> $Project.key
#echo "QMMM-ELECTROSTATICS ESPF" >> ${Project}_nolink.key
echo "QMMM-microiteration ON" >> ${Project}_nolink.key

cat chargemm >> $Project.key
cat chargemm >> ${Project}_nolink.key
cat qm >> mm
for i in $(eval echo "{1..$(($numchr-$xx-$fx))}"); do
   active=`head -n $i mm | tail -n1 | awk '{ print $2 }'`
   echo "ACTIVE $active" >> $Project.key
   echo "ACTIVE $active" >> ${Project}_nolink.key
done
echo "ACTIVE $(($final+1))" >> $Project.key
rm -f qm mm col

#
# Adding the link atom to the xyz and key files using tinker
#
mv ${Project}_nolink.key $Project-tk.key
$tinkerdir/xyzedit $Project-tk.xyz << EOF
20
0
EOF

cp $Project-tk.xyz_2 ${Project}_ESPF.xyz
cp $Project.key ${Project}_ESPF.key

multchain=`grep "MultChain" ../Infos.dat | awk '{ print $2 }'`
if [[ $multchain == "YES" ]]; then
   answer="b"
   while [[ $answer != "y" && $answer != "n" ]]; do
      echo ""
      echo "*****************************************************"
      echo " This seems to be a multi-chain protein. So, if this"
      echo " is a dimer calculation (one chromophore per chain)"
      echo " you must provide later the full path of the charges file"
      echo ""
      echo " Is this a dimer calculation? (y/n)"
      echo "*****************************************************"
      read answer
   done
   if [[ $answer == "y" ]]; then
      ../update_infos.sh "Dimer" "YES" ../Infos.dat
   else
      ../update_infos.sh "Dimer" "NO" ../Infos.dat
   fi
else
   answer="n"
   ../update_infos.sh "Dimer" "NO" ../Infos.dat
fi

#
# Single point molcas calculation to fit the charges of the chromophore
# (self explaining echos)
#
if [[ $answer == "n" ]]; then

   echo " **********************************************************************"
   echo ""
   echo " The charges of the Chromophore will be fit now by using the ESPF model"
   echo ""
   echo " **********************************************************************"

   cp $templatedir/molcas.slurm.sh molcas-job.sh
   NOME=${Project}_ESPF
   sed -i "s|NOMEPROGETTO|$NOME|" molcas-job.sh
   no=$PWD
   sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
   sed -i "s|MEMTOT|23000|" molcas-job.sh
   sed -i "s|MEMORIA|20000|" molcas-job.sh
   sed -i "s|hh:00:00|100:00:00|" molcas-job.sh

   cp $templatedir/ASEC/templateSP ${Project}_ESPF.input

   state=10
   while  [[ $state -ne 0 && $state -ne 1 && $state -ne 2 && $state -ne 3 && $state -ne 4 && $state -ne 5 && $state -ne 6 && $state -ne 7 ]]; do
      echo ""
      echo " What is the electronic state to be optimized?" 
      echo " (1 corresponds to the grouns state)"
      echo ""
      echo ""
      read state
   done

   if [[ $state -gt 1 ]]; then
      Ras2=
      sed -i "s/\*   ciroot=3 3 1/\ \ \ \ ciroot=$state $state 1/" ${Project}_ESPF.input
      sed -i "/ciroot=/a\ \ \ \ rlxroot=$state" ${Project}_ESPF.input
   fi
   ../update_infos.sh "Electronic_state" $state ../Infos.dat


   option=0
   while [[ $option -ne 1 && $option -ne 2 ]]; do
      echo ""
      echo " Define the active space. Please select one option:"
      echo ""
      echo " 1) Default option: Neutral flavin with a CASSCF comprising 10 electrons in 10 orbitals."
      echo " 2) Customize the active space"
      echo ""
      echo ""
      read option
   done

   if [[ $option -eq 1 ]]; then
      ../update_infos.sh "Active_space" "default" ../Infos.dat
   else
      ../update_infos.sh "Active_space" "custom" ../Infos.dat
      echo ""
      echo " How many active electrons:"
      echo ""
      read actelec
      echo ""
      echo " How many active orbitals:"
      echo ""
      read actorb
      echo ""
      echo " How many inactive orbitals:"
      echo ""
      read inactorb
      sed -i "s|nActEl=10 0 0|nActEl=$actelec 0 0|" ${Project}_ESPF.input
      sed -i "s|Ras2=10|Ras2=$actorb|" ${Project}_ESPF.input
      sed -i "s|Inactive=62|Inactive=$inactorb|" ${Project}_ESPF.input

      ../update_infos.sh "Active_electrons" "$actelec" ../Infos.dat
      ../update_infos.sh "Active_orbitals" "$actorb" ../Infos.dat
      ../update_infos.sh "Inactive_orbitals" "$inactorb" ../Infos.dat
   fi
   sed -i "s|PARAMETRI|${prm}|" ${Project}_ESPF.input
   sed -i "s|VDZ|VDZP|" ${Project}_ESPF.input

   echo ""
   echo ""
   echo " ***********************************************************************"
   echo ""
   echo " A single point calculation submitted (CASSCF/ANO-L-VDZP)"
   echo " to generate the ESPF charges."
   echo ""
   echo " ***********************************************************************"
   echo ""
   echo ""
else
#
# This part needs to be re-done!!
#
   filecheck=b
   while [[ $filecheck != "OK" ]]; do
      echo ""
      echo "**********************************************************"
      echo " Please provide the full path of the charges file,"
      echo " i.e. /home/.../new_charges"
      echo "**********************************************************"
      echo ""
      read chargepath
      if [[ -f $chargepath ]]; then
         filecheck="OK"
      fi
   done

   if [[ $option -eq 1 ]]; then
      ../update_infos.sh "iLOV_RESP" "iLOV" ../Infos.dat
   else
      ../update_infos.sh "iLOV_RESP" "iLOV_chain" ../Infos.dat
   fi

   cp $chargepath new_charges
   num=$(($numchromo/2))
   head -n $num new_charges > tempo
   head -n $num new_charges >> tempo
   echo "" >> tempo
   mv tempo new_charges
   cp new_charges ../
fi

sbatch molcas-job.sh

#
# Message to the user
#
cd ..

cp $templatedir/ASEC/fitting_ESPF.sh .
./update_infos.sh "Next_script" "fitting_ESPF.sh" Infos.dat

echo ""
echo "**********************************************************"
echo " Wait for the Molcas single point calculation to be done,"
echo " then execute: fitting_ESPF.sh"
echo "**********************************************************"
echo ""

