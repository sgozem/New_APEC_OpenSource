#!/bin/bash
#
# Reading information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
multichain=`grep "MultChain" Infos.dat | awk '{ print $2 }'`
step=`grep "Step" Infos.dat | awk '{ print $2 }'`
charge=`grep "Init_Charge" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
dimer=`grep "Dimer" Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`

#
# Having parametrized the charges of the chromophore using the ESPF model
# we can go back to the Minimize_${Project} folder to continue with the 
# energy minimization
#
# Retriving the pdb of the chromophore and rtp file
#
cd Minimize_${Project}
rm -rf $amber.ff *.itp *.top $Project.gro
cp ../Chromophore/${chromophore}.pdb .
cp -r $templatedir/$amber.ff .
cd $amber.ff/


option=0
while [[ $option -ne 1 && $option -ne 2 && $option -ne 3 ]]; do
   echo ""
   echo " Please select the Flavin model to use:"
   echo ""
   echo " 1) Quinone"
   echo " 2) Semi-quinone"
   echo " 3) Hydro-quinone"
   echo ""
   read option
done
../../update_infos.sh "Redox" $option ../../Infos.dat

fmnfad="NONE"
while [[ $fmnfad != "FMN" && $fmnfad != "FAD" ]]; do
   echo ""
   echo " Please select if the tail corresponds to FMN or FAD (just type FMN or FAD)"
   echo ""
   echo ""
   read fmnfad
done
../../update_infos.sh "Tail" $fmnfad ../../Infos.dat


if [[ $option -eq 1 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/manchester_FMN_rtp new_rtp
#else
#   echo ""
#   echo "*******************************************************"
#   echo ""
#   echo " Implement this parametrizations from Manchester: "
#   echo " http://research.bmh.manchester.ac.uk/bryce/amber/"
#   echo ""
#   echo "*******************************************************"
#   exit 0
fi
if [[ $option -eq 2 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/manchester_FMNH_rtp new_rtp
fi
if [[ $option -eq 3 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/manchester_FMNH2_rtp new_rtp
fi

#cp ../../ESPF_charges/new_rtp .
cat new_rtp >> aminoacids.rtp
cd ..

#
# Modifying the pdb file to put the chromophore between the protein and the waters
#
wat=`grep "DOWSER_wat" ../Infos.dat | awk '{ print $2 }'`
nwat=$(($wat+$wat+$wat))
if [[ $nwat -gt 0 ]]; then
   numchromo=`grep -c " CHR " ${chromophore}.pdb`
   lineas=`wc -l $Project.pdb | awk '{ print $1 }'`
   head -n $(($lineas-$nwat)) $Project.pdb > proteina
   tail -n $nwat $Project.pdb > tempwat
   cat ${chromophore}.pdb >> proteina
   head -n $(($lineas+$numchromo-$nwat)) proteina > proteina2
   cat tempwat >> proteina2
   mv proteina2 $Project.pdb
   rm proteina
else
   cat ${chromophore}.pdb >> $Project.pdb
fi

#
# pdb2gmx is the Gromacs utility for generating gro files and topologies
#
$gropath/gmx pdb2gmx -f $Project.pdb -o $Project.gro -p $Project.top -ff $amber -water tip3p 2> grolog
checkgro=`grep 'Writing coordinate file...' grolog`
   if [[ -z $checkgro ]]; then
      echo " An error occurred during the execution of pdb2gmx. Please look into grolog file"
      echo " No further operation performed. Aborting..."
      echo ""
      exit 0
   else
      echo " new.gro and its topology were successfully generated"
      echo ""
      rm grolog
   fi

#
# Generating the ndx file for moving sidechains
#
echo '2' > choices.txt
echo 'q' >> choices.txt
$gropath/gmx make_ndx -f $Project.gro -o $Project.ndx < choices.txt
grep -n 'OW' $Project.gro | cut -d : -f 1 > temporal
awk '$1=$1-2' temporal > oxywat
echo "[ Group1 ]" >> $Project.ndx
cat oxywat >> $Project.ndx
rm oxywat choices.txt

backb=0
#while  [[ $backb -ne 0 && $backb -ne 1 ]]; do
       echo " *******************************************************"
       echo ""
       echo " An energy minimization of the protein side-chains and"
       echo " hydrogens of the chromophore will be performed before"
       echo " embedding in the Solvent Box. The backbone will not be"
       echo " relaxed for the moment."
       echo ""
#       echo " Please type 1 if you want to relax the backbone, 0 otherwise"
       echo " *******************************************************"
       echo ""
#       read backb
#done
sleep 10

cp $templatedir/ASEC/ndx-maker_mod.sh .
./ndx-maker_mod.sh $Project 5 $backb

#
# Runing the MM side chains energy minimization in the loging node.
# In order to do this we divided the minimization steps in batches
# of 1000 steps. Otherwise it would be killed bu the system
#
conver=10
iter=1
while [[ conver -ne 0 ]]; do
   $gropath/gmx grompp -f standard-EM.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr -maxwarn 2

   echo ""
   echo " Please wait, minimizing, batch $iter of 1000 steps"
   echo ""

   $gropath/gmx mdrun -s $Project.tpr -o $Project.trr -x $Project.xtc -c final-$Project.gro 2> grolog

   echo '0' > choices.txt
   echo 'q' >> choices.txt
   $gropath/gmx trjconv -pbc nojump -f final-$Project.gro -s $Project.gro -o ${Project}_No_jump.gro < choices.txt


   if grep -q "Steepest Descents did not converge to Fmax" md.log; then
      mkdir Iter_$iter
      mv ener.edr $Project.gro final-$Project.gro $Project.tpr $Project.trr md.log mdout.mdp Iter_$iter
      cp ${Project}_No_jump.gro Iter_$iter
      mv ${Project}_No_jump.gro $Project.gro
      iter=$(($iter+1))
   else
      if grep -q "Steepest Descents converged to Fmax" md.log; then
         conver=0
         echo ""
         echo " MM energy minimization seems to finish properly."
         echo ""
         mv final-$Project.gro backup_final-$Project.gro
         cp ${Project}_No_jump.gro final-$Project.gro
      else
         echo ""
         echo " There is a problem with the energy minimization. Please check it."
         echo ""
         ans="b"
         while [[ $ans != "y" && $ans != "n" ]]; do
            echo "******************************************************************"
            echo ""
            echo " Do you still want to continue? (y/n)"
            echo ""
            echo "******************************************************************"
            read ans
         done
         if [[ $ans == "n" ]]; then
            exit 0
         else
            conver=0
            mv final-$Project.gro backup_final-$Project.gro
            cp ${Project}_No_jump.gro final-$Project.gro
         fi
      fi
   fi
done

cd ..

#
# When the minimization of the side chains is done the whole system will be
# embbeded in a solvent box and the whole system will be minimized. 
#
#continua="a"
#while [[ $continua != "y" && $continua != "n" ]]; do
echo " ****************************************************************"
echo ""
echo " The protein will be embbeded in a solvent box. Then it will be" 
echo " energetically minimized relaxing also the backbone." 
#echo " Continue? (y/n)"
echo ""
echo " ****************************************************************"
#read continua
#done
sleep 10

#if [[ $continua == "n" ]]; then
#   echo ""
#   echo " Aborting..."
#   echo ""
#   exit 0
#fi

echo " ****************************************************************"
echo ""
echo " Please define the size of the cubic box in nanometers (i.e 7.0)."
echo " For a single proteins 7.0 is normally ok."
echo " For a dimer, maybe 10.0 is fine"
echo ""
echo " ****************************************************************"
read box

relaxpr="y"
./update_infos.sh "Relax_protein" "$relaxpr" Infos.dat

if [[ $relaxpr == "y" ]]; then
   relaxbb="y"
#   while [[ $relaxbb != y && $relaxbb != n ]]; do
#      echo ""
#      echo " Relax backbone? (y/n)"
#      echo ""
#      read relaxbb
#   done
   ./update_infos.sh "Relax_backbone" "$relaxbb" Infos.dat
else
   ./update_infos.sh "Relax_backbone" "n" Infos.dat 
fi

mkdir Dynamic
cd Dynamic
mkdir Box
cd ..
cd Minimize_$Project
cp -r $amber.ff *.itp $Project.top residuetypes.dat ../Dynamic/Box
cp final-$Project.gro ../Dynamic/Box/$Project.gro
cd ..
cd Dynamic/Box

$gropath/gmx editconf -f $Project.gro -bt cubic -box $box $box $box -o ${Project}_box_init.gro -c
$gropath/gmx solvate -cp ${Project}_box_init.gro -cs spc216.gro -o ${Project}_box_sol_init.gro -p $Project.top >& genbox.out

#
# Excluding the water molecules from the solvent box that can be added
# very close to the protein (like in cavities). Otherwise, it can bring 
# problems in the energy minimization).
# The water molecules added to fill the box are labeled as SOL, while the ones coming from dowser are HOH.
# This avoids the possibility of removing internal waters originally coming from the pdb.
#
cp ${Project}_box_sol_init.gro Interm.gro

selection="((same residue as all within 0.5 of protein) and resname SOL) or ((same residue as all within 2 of resname CHR) and resname SOL)"

# TCL script for VMD: open file, apply selection, save the serial numbers into a file
#
echo -e "mol new Interm.gro" > removewat.tcl
echo -e "mol delrep 0 top" >> removewat.tcl
line1="set sele [ atomselect top \"$selection\" ]" 
echo -e "$line1" >> removewat.tcl
echo -e 'set numbers [$sele get serial]' >> removewat.tcl
line2="set filename Interwat"
echo -e "$line2" >> removewat.tcl
echo -e 'set fileId [open $filename "w"]' >> removewat.tcl
echo -e 'puts -nonewline $fileId $numbers' >> removewat.tcl
echo -e 'close $fileId' >> removewat.tcl
echo -e "exit" >> removewat.tcl
vmd -e removewat.tcl -dispdev text
rm removewat.tcl
echo ""
echo ""
echo " Please wait ..."

#
# Excluding those waters from the gro file
#
cp Interm.gro Interm_yoe.gro

numinit=`head -n2 $Project.gro | tail -n1 | awk '{ print $1 }'`
col=`awk '{print NF}' Interwat`
cont=0
for i in $(eval echo "{1..$col}")
do
   rem=`expr $i % 3`
   indx=`awk -v j=$i '{ print $j }' Interwat`
   if [[ $indx -gt $numinit ]]; then
      cont=$(($cont+1))
      if [[ $rem -eq 1 ]]; then
         sed -i "/  OW$indx /d" Interm.gro
         sed -i "/  OW $indx /d" Interm.gro
         sed -i "/  OW  $indx /d" Interm.gro
         sed -i "/  OW   $indx /d" Interm.gro
         sed -i "/  OW    $indx /d" Interm.gro
      else
         if [[ $rem -eq 2 ]]; then
            sed -i "/ HW1$indx /d" Interm.gro
            sed -i "/ HW1 $indx /d" Interm.gro
            sed -i "/ HW1  $indx /d" Interm.gro
            sed -i "/ HW1   $indx /d" Interm.gro
            sed -i "/ HW1    $indx /d" Interm.gro
         else
            sed -i "/ HW2$indx /d" Interm.gro
            sed -i "/ HW2 $indx /d" Interm.gro
            sed -i "/ HW2  $indx /d" Interm.gro
            sed -i "/ HW2   $indx /d" Interm.gro
            sed -i "/ HW2    $indx /d" Interm.gro
         fi
      fi
   fi
done

numatm=`head -n2 Interm.gro | tail -n1 | awk '{ print $1 }'`
newnum=$(($numatm-$cont))
head -n1 Interm.gro > ${Project}_box_sol.gro
echo "$newnum" >> ${Project}_box_sol.gro
tail -n$(($numatm-$cont+1)) Interm.gro >> ${Project}_box_sol.gro
#rm Interm.gro

#
# generating the ndx and top files for the energy minimization of the whole system
# (solvent box plus protein)
#
echo '2' > choices.txt
echo 'q' >> choices.txt
$gropath/gmx make_ndx -f ${Project}_box_sol.gro -o ${Project}_box_sol.ndx < choices.txt
rm choices.txt

addwat=`tail -n1 $Project.top | awk '{ print $2 }'`
wattop=$((($addwat*3-$cont)/3))

lines=`wc -l $Project.top | awk '{ print $1 }'`
cp $Project.top tempo
head -n$(($lines-1)) tempo > ${Project}_box_sol.top
echo -e "SOL              $wattop" >> ${Project}_box_sol.top
rm tempo

mkdir ../Minimization
cp $templatedir/ASEC/min_sol.mdp ../Minimization
cp -r $amber.ff *.itp residuetypes.dat ../Minimization
cp ${Project}_box_sol.gro ${Project}_box_sol.ndx ${Project}_box_sol.top ../Minimization

cd ../Minimization

#
# Defining the "GroupDyna" group of the ndx file to be fixed during the 
# MM energy minimization
#
chratoms=`head -n1 ../../Chromophore/$chromophore.xyz | awk '{ print $1 }'`

if [[ $relaxpr == y ]]; then
   if [[ $relaxbb == n ]]; then
      selection="backbone or (resname CHR and not hydrogen)"
      # TCL script for VMD: open file, apply selection, save the serial numbers into a file
      #
      echo -e "mol new ${Project}_box_sol.gro type gro" > ndxsel.tcl
      echo -e "mol delrep 0 top" >> ndxsel.tcl
      riga1="set babbeo [ atomselect top \"$selection\" ]"
      echo -e "$riga1" >> ndxsel.tcl
      echo -e 'set noah [$babbeo get serial]' >> ndxsel.tcl
      riga3="set filename dinabb"
      echo -e "$riga3" >> ndxsel.tcl
      echo -e 'set fileId [open $filename "w"]' >> ndxsel.tcl
      echo -e 'puts -nonewline $fileId $noah' >> ndxsel.tcl
      echo -e 'close $fileId' >> ndxsel.tcl
      echo -e "exit" >> ndxsel.tcl
      vmd -e ndxsel.tcl -dispdev text

      num=`awk '{print NF}' dinabb`
      for i in $(eval echo "{1..$num}")
      do
        awk -v j=$i '{ print $j }' dinabb >> dina
      done
   else
#
# Selecting the atoms of the chromophore plus the fixed and LQ, LM
# to be fixed during the MD
#
      grep -n 'CHR ' ${Project}_box_sol.gro | cut -d : -f 1 > temporal1
      awk '$1=$1-2' temporal1 > dina
      rm -f dina2
      for i in $(eval echo "{1..$chratoms}"); do
         atmtype=`head -n $(($i+2)) ../../Chromophore/$chromophore.xyz | tail -n1 | awk '{ print $6 }'`
         if [[ $atmtype == "QM" || $atmtype == "MM" || $atmtype == "LM" || $atmtype == "LQ" ]]; then
            head -n $i dina | tail -n1 | awk '{ print $1 }' >> dina2
         fi
      done
      mv dina2 dina
      num=`grep -c "CHR " ${Project}_box_sol.gro | awk '{ print $1 }'`
   fi
   if [[ $dimer == "YES" ]]; then
      chrnum=`grep -c "CHR " ${Project}_box_sol.gro | awk '{ print $1 }'`
      head -n $(($num-$chrnum+$chrnum/2)) dina > dinadimer
      mv dinadimer dina 
   fi
fi

tr '\n' ' ' < dina > dyna
echo ":set tw=75" > shiftline.vim
echo ":normal gqw" >> shiftline.vim
echo ":x" >> shiftline.vim
vim -es dyna < shiftline.vim
echo "[ GroupDyna ]" >> ${Project}_box_sol.ndx
cat ${Project}_box_sol.ndx dyna > last.ndx
mv last.ndx ${Project}_box_sol.ndx

sed -i "s/;freezegrps = GroupDyna/freezegrps = GroupDyna/g" min_sol.mdp
sed -i "s/;freezedim = Y Y Y/freezedim = Y Y Y/g" min_sol.mdp

#
# Generating the tpr for the MM energy minimization
#
$gropath/gmx grompp -f min_sol.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr -maxwarn 2

#   
# Adding ions to the Solvent Box for neutralizing the total charge of the system
# (see self explaning echoes)
#

pairs=-1
while [[ $pairs -lt 0 ]]; do

   if [[ $charge -ne 0 ]]; then
      echo ""
      echo " *************************************************************"
      echo "  The total charge of the system is ${charge}, which will be"
      echo "  neutralized by adding Na or Cl ions."
      echo "  But, if you want to add extra pair of ions to the system"
      echo "  to mimmic the experimental conditions, please specify"
      echo "  how many pairs of ions (NA and CL) to add." 
      echo "  Type \"0\" otherwise."
      echo ""
      read pairs
   fi

   if [[ $charge -eq 0 ]]; then
      echo ""
      echo " *************************************************************"
      echo "  The total charge of the system is ZERO, no Na or Cl ions"
      echo "  will be added to neutralize the system."
      echo "  But, if you want to add extra pair of ions to the system"
      echo "  to mimmic the experimental conditions, please specify"
      echo "  how many pairs of ions (NA and CL) to add." 
      echo "  Type \"0\" otherwise."
      echo ""
      read pairs
   fi
done

if [[ $charge -ne 0 || $pairs -ne 0 ]]; then

   mkdir Add_Ion
   mv ${Project}_box_sol.tpr ${Project}_box_sol.ndx ${Project}_box_sol.top ${Project}_box_sol.gro Add_Ion
   cd Add_Ion
   numpro=`head -n2 ../../Box/$Project.gro | tail -n1 | awk '{ print $1 }'`

#
# This while is used to ensure that the ions will not be added inside the proteins
#
   res=0
   seed=111
   count=0
   while [[ $res -eq 0 ]]; do
      replace=0
      replacectl=0
      if [[ $charge -lt 0 ]]; then
         pcharge=$(echo "-1*$charge" | bc)
         pcharge=$(($pcharge+$pairs)) 
         if [[ -f back_${Project}_box_sol.top ]]; then
            cp back_${Project}_box_sol.top ${Project}_box_sol.top
         else
            cp ${Project}_box_sol.top back_${Project}_box_sol.top
         fi
         $gropath/gmx genion -seed $seed -s ${Project}_box_sol.tpr -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -pname NA -pq 1 -np $pcharge -nname CL -nq -1 -nn $pairs -o ${Project}_box_sol_ion.gro 2> addedions << EOF
15
EOF
         
#         ../../../update_infos.sh "Added_Ions" "${pcharge}_NA" ../../../Infos.dat
         ../../../update_infos.sh "Added_NAs" "$pcharge" ../../../Infos.dat
	 ../../../update_infos.sh "Added_CLs" "$pairs" ../../../Infos.dat

         lin=`grep -c "Replacing solvent molecule" addedions`
         for i in $(eval echo "{1..$lin}")
         do
            replace=`grep "Replacing solvent molecule" addedions | head -n $i | tail -n1 | awk '{ print $6 }' | sed 's/[^0-9]//g'`
         if [[ $replace -le $numpro ]]; then
            replacectl=$(($replacectl+1))
         fi
         done
      fi
      if [[ $charge -ge 0 ]]; then
         if [[ -f back_${Project}_box_sol.top ]]; then
            cp back_${Project}_box_sol.top ${Project}_box_sol.top
         else
            cp ${Project}_box_sol.top back_${Project}_box_sol.top
         fi
         charge=$(($charge+$pairs))
         $gropath/gmx genion -seed $seed -s ${Project}_box_sol.tpr -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -nname CL -nq -1 -nn $charge -pname NA -pq 1 -np $pairs -o ${Project}_box_sol_ion.gro 2> addedions << EOF
15
EOF
#         ../../../update_infos.sh "Added_Ions" "${charge}_CL" ../../../Infos.dat
         ../../../update_infos.sh "Added_CLs" "$charge" ../../../Infos.dat
         ../../../update_infos.sh "Added_NAs" "$pairs" ../../../Infos.dat

         lin=`grep -c "Replacing solvent molecule" addedions`
         for i in $(eval echo "{1..$lin}")
         do
            replace=`grep "Replacing solvent molecule" addedions | head -n $i | tail -n1 | awk '{ print $6 }' | sed 's/[^0-9]//g'`
         if [[ $replace -le $numpro ]]; then
            replacectl=$(($replacectl+1))
         fi
         done
   
      fi
#
#  This while is to ensure that the SOL group has been selected by "genion" code for placing the ions.
#  In other systems it could be defferent from 13 and it should be corrected.
#
      soll=b
      while [[ $soll != "y" && $soll != "n" ]]; do
         echo ""
         echo " Was selected the \"SOL\" group for adding the ions? (y/n)"
         echo ""
         read soll
         if [[ $soll == "n" ]]; then
            echo ""
            echo " Modify this script in order to sellect the right"
            echo " number of the \"SOL\" group"
            echo " Terminating ..."
            echo ""
            exit 0
         fi
      done

      if [[ $replacectl -eq 0 ]]; then
         res=10

         cp ${Project}_box_sol_ion.gro ../${Project}_box_sol.gro
         cp ${Project}_box_sol.top ../${Project}_box_sol.top

         cd ..
         echo '2' > choices.txt
         echo 'q' >> choices.txt
         $gropath/gmx make_ndx -f ${Project}_box_sol.gro -o ${Project}_box_sol.ndx < choices.txt
         rm choices.txt
         echo "[ GroupDyna ]" >> ${Project}_box_sol.ndx
         cat ${Project}_box_sol.ndx dyna > last.ndx
         mv last.ndx ${Project}_box_sol.ndx
         rm dyna
      else
         seed=$(($seed+3))
      fi
   done
fi

if [[ $charge -eq 0 && $pairs -eq 0 ]]; then
   ../../update_infos.sh "Added_CLs" 0 ../../Infos.dat
   ../../update_infos.sh "Added_NAs" 0 ../../Infos.dat
fi

#
# Runing the MM energy minimization. Again it will run in the loging node dividing the number of steps
# in batches of 1000 steps to avoid being killed the the system
# (See the self explaning echoes)
#
conver=10
iter=1
if [[ $dimer == "YES" ]]; then
   sed -i "s/emtol                   = 100/emtol                   = 500/" min_sol.mdp
fi
while [[ conver -ne 0 ]]; do
   $gropath/gmx grompp -f min_sol.mdp -c ${Project}_box_sol.gro -n ${Project}_box_sol.ndx -p ${Project}_box_sol.top -o ${Project}_box_sol.tpr -maxwarn 2

   echo ""
   echo " Please wait, minimizing, batch $iter of 1000 steps"
   echo ""

   $gropath/gmx mdrun -s ${Project}_box_sol.tpr -o ${Project}_box_sol.trr -x ${Project}_box_sol.xtc -c final-${Project}_box_sol.gro 2> grolog

   if grep -q "Steepest Descents did not converge to Fmax" md.log; then
      mkdir Iter_$iter
      mv ener.edr ${Project}_box_sol.gro ${Project}_box_sol.tpr ${Project}_box_sol.trr md.log mdout.mdp Iter_$iter
      cp final-${Project}_box_sol.gro Iter_$iter
      mv final-${Project}_box_sol.gro ${Project}_box_sol.gro
      iter=$(($iter+1))
   else 
      if grep -q "Steepest Descents converged to Fmax" md.log; then
         conver=0  
         echo ""
         echo " MM energy minimization seems to finished properly."
         echo ""
         cd ../../
#         cp $templatedir/ASEC/ESPF_charges.sh .
#         ./update_infos.sh "Next_script" "ESPF_charges.sh" Infos.dat
         cp $templatedir/ASEC/MD_NPT.sh .
         ./update_infos.sh "Next_script" "MD_NPT.sh" Infos.dat
         ./update_infos.sh "MD_ensemble" "NVT" Infos.dat
         echo ""
         echo "**********************************************************"
         echo ""
#         echo " Pleasde execute now: ESPF_charges.sh"
         echo " Pleasde execute now: MD_NPT.sh"
         echo ""
         echo "**********************************************************"
         echo ""
      else
         echo ""
         echo "*********************************************************************************"
         echo ""
         echo " There is a problem with the energy minimization. Please check it. Terminating..."
         echo ""
         echo "*********************************************************************************"
         exit 0
      fi
   fi
done

