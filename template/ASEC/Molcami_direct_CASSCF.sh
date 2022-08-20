#!/bin/bash
#
# Retrieving the needed information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
uphess=`grep "Update_hessian" Infos.dat | awk '{ print $2 }'`
solvent=`grep "SolventBox" Infos.dat | awk '{ print $2 }'`
shell=`grep "Shell" Infos.dat | awk '{ print $2 }'`
stationary=`grep "Stationary" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
state=`grep "Electronic_state" Infos.dat | awk '{ print $2 }'`
ilov=`grep "iLOV_RESP" Infos.dat | awk '{ print $2 }'`
chargechr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`

#
# Copying and running the script for Tinker files preparation
#
#
# Creating the directory calculations and copying all the needed files
# If calculations already exists, the script exits with a message
#
echo ""
echo " Using $Project.xyz and $prm.prm"
echo ""
if [ -d calculations ]; then
   echo " Folder calculations already exists. Please check it and remove if necessary"
   echo " Terminating..."
   echo ""
   exit 0
fi

#
# Creating key file for the ASEC QMMM calculations
#
echo "parameters $prm" > $Project.key
echo "parameters $prm" > ${Project}_nolink.key
echo "" >> $Project.key
echo "" >> ${Project}_nolink.key
chratoms=`head -n1 Chromophore/$chromophore.xyz | awk '{ print $1 }'`
xx=0
fx=0
final=`head -n1 $Project-tk.xyz | awk '{ print $1 }'`
init=$(($final-$chratoms))
rm -f mm
rm -f qm
rm -f chargemm
rm -f chargexx
rm -f chargefx
rm -f chargefxx
touch mm
touch qm
for i in $(eval echo "{1..$chratoms}"); do
   num=`head -n $(($init+1+$i)) $Project-tk.xyz | tail -n1 | awk '{ print $1 }'`
   charge=`head -n $i new_charges | tail -n1 | awk '{ print $1 }'`
   atmtype=`head -n $(($i+2)) Chromophore/$chromophore.xyz | tail -n1 | awk '{ print $6 }'`
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
      echo " When FAD chromophore is computed, like in this case, just the LM atom"
      echo " will be considered in the QMMM calculations as MM atoms. The rest of the"
      echo " tail will be considered as ASEC points."
      echo ""
      echo " Do you agree? (y/n)"
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

#
# The charges of the remaining MM atoms are rounded to $carga to ensure that the
# total charge of the QM+MM part is the desired charge.
#
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

echo "QMMM $(($chratoms-$xx-$fx+1))" >> $Project.key
echo "QMMM $(($chratoms-$xx-$fx))" >> ${Project}_nolink.key

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
for i in $(eval echo "{1..$(($chratoms-$xx-$fx))}"); do
   active=`head -n $i mm | tail -n1 | awk '{ print $2 }'`
   echo "ACTIVE $active" >> $Project.key
   echo "ACTIVE $active" >> ${Project}_nolink.key
done
echo "ACTIVE $(($final+1))" >> $Project.key
rm -f qm mm col
#chargemm chargemm0

#
# Adding the link atom to the tinker xyz and key files
#
mv ${Project}_nolink.key $Project-tk.key
$tinkerdir/xyzedit $Project-tk.xyz << EOF
20
0
EOF

if [[ -f $Project-tk.xyz_2 ]] ; then
   if grep -q "HLA " $Project-tk.xyz_2; then
      if [[ $Molcami_OptSCF != "YES" ]]; then
#         ./update_infos.sh "Shell" $(($shell+1)) Infos.dat
         ./update_infos.sh "Molcami_OptSCF" "YES" Infos.dat
      fi
#
# Reformating the final xyz tinket file
#
cat > reformat.f << YOE
      Program reformat
      implicit real*8 (a-h,o-z)
      character line3*3,line74*74
      open(1,file='$Project-tk.xyz_2',status='old')
      open(2,file='$Project-tk.xyz_new',status='unknown')
      read(1,*)
      write(2,'(i6)')$(($final+1))
      do i=1,$(($final+1))
           read(1,'(i5,1x,A,1x,A)')nume,line3,line74
           write(2,'(i6,2x,A,1x,A)')nume,line3,line74
      enddo
      end
YOE
      gfortran reformat.f -o reformat.x
      ./reformat.x
      rm $Project-tk.xyz_2
      mv $Project-tk.xyz_new $Project-tk.xyz_2
      rm $Project-tk.key
#      rm $Project-tk.input
   else
      echo "*********************************"
      echo " It seems to be something wrong"
      echo " adding the Link atom."
      echo " Aborting ..."
      echo "*********************************"
      exit 0
   fi
else
   echo "*********************************"
   echo " It seems to be something wrong"
   echo " adding the Link atom."
   echo " Aborting ..."
   echo "*********************************"
   exit 0
fi

#
# In the calculations folder will be run all the QM/MM geometries optimizations
# and CASPT2 calculations
#
mkdir calculations
newdir=${Project}_VDZP_Opt
mkdir calculations/${newdir}

cp $Project.key calculations/${newdir}/${newdir}.key

#
# This section is for inserting the optimized cromophore structure of the previous 
# Step into the new -tk.xyz in order to avoid the rounding arising from the tk to gro
# conversion and the new position of the Link atom.
#
if [[ $Step -gt 0 ]]; then
   cp ../Step_$(($Step-1))/calculations/${newdir}/${newdir}.Final.xyz ${newdir}.Final_last.xyz

   cp $templatedir/ASEC/Ins_tk_tk_Molcami.f .
   cp $Project-tk.xyz_2 new_tk.xyz
   cp ${newdir}.Final_last.xyz last_tk.xyz

   sed -i "s|numero|$(($shell+1))|g" Ins_tk_tk_Molcami.f

   grep "QM \|MM \|LA " $Project.key | grep -v "QMMM" | awk '{ print $2 }' > qmmmatoms
   sort -n -o sorted qmmmatoms
   mv sorted qmmmatoms

   fin=`wc -l qmmmatoms | awk '{ print $1 }'`
   sed -i "s/finall/$fin/g" Ins_tk_tk_Molcami.f

   gfortran Ins_tk_tk_Molcami.f -o Ins_tk_tk_Molcami.x
   ./Ins_tk_tk_Molcami.x
   mv Insert-tk.xyz calculations/${newdir}/${newdir}.xyz
else
   mv $Project-tk.xyz_2 calculations/${newdir}/${newdir}.xyz
fi
calcul=1

#
# Collecting data from the previos step
#
if [[ -a ${newdir}.JobIph_old ]]; then
   cp ${newdir}.JobIph_old calculations/${newdir}/${newdir}.JobIph
fi
if [[ -a ${newdir}.RasOrb ]]; then
   cp ${newdir}.RasOrb calculations/${newdir}/${newdir}.RasOrb
fi
if [[ -a ${newdir}.Hessian_old ]]; then
   cp ${newdir}.Hessian_old calculations/${newdir}/${newdir}.Hessian_old
fi
if [[ -a ${newdir}.Hessian_QMMM ]]; then
   cp ${newdir}.Hessian_QMMM calculations/${newdir}/${newdir}.Hessian_old
fi

cp $templatedir/molcas.slurm.sh calculations/${newdir}/molcas-job.sh


if [[ $retstereo == "nAT" ]]; then
   cp $templatedir/templateOPTSCFneu calculations/${newdir}/templateOPTSCF 
else
 if [[ $calcul == 1 ]]; then
   cp $templatedir/ASEC/template_CASSCF_min calculations/${newdir}/template_CASSCF
   sed -i "s|> COPY \$InpDir/\$Project.Espf.Data \$WorkDir|*> COPY \$InpDir/\$Project.Espf.Data \$WorkDir|" calculations/${newdir}/template_CASSCF
   if [[ $state -gt 1 ]]; then
      sed -i "/Ras2=/a\ \ \ \ ciroot=$state $state 1" calculations/${newdir}/template_CASSCF
      sed -i "/ciroot=/a\ \ \ \ rlxroot=$state" calculations/${newdir}/template_CASSCF
      if [[ $Step -eq 1 ]]; then
         sed -i "/\ \ \ \ cirestart/d" calculations/${newdir}/template_CASSCF
      fi
   fi
 else
   cp $templatedir/ASEC/template_CASSCF_TS calculations/${newdir}/template_CASSCF
   if [[ $Step -eq 0 ]]; then
      sed -i "/JobIph_old/d" calculations/${newdir}/template_CASSCF
      sed -i "/\ \ \ \ JobIph/d" calculations/${newdir}/template_CASSCF
      sed -i "/\ \ \ \ cirestart/d" calculations/${newdir}/template_CASSCF
   fi
 fi
fi

#
# Submiting the QM/MM geometry optimization
#
cd calculations
cd ${newdir}/

NOME=${Project}_VDZP_Opt
sed -i "s|NOMEPROGETTO|$NOME|" molcas-job.sh

no=$PWD
sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|23000|" molcas-job.sh
sed -i "s|MEMORIA|20000|" molcas-job.sh
sed -i "s|hh:00:00|120:00:00|" molcas-job.sh

sed -i "s|PARAMETRI|${prm}|" template_CASSCF

if [[ $uphess == NO ]]; then
   sed -i "/RUNOLD/d" template_CASSCF
   sed -i "/OLDF/d" template_CASSCF
fi

mv template_CASSCF ${newdir}.input 
sed -i "s|VDZ|VDZP|" ${newdir}.input

actspace=`grep "Active_space" ../../Infos.dat | awk '{ print $2 }'`
if [[ $actspace == "custom" ]]; then
   actelec=`grep "Active_electrons" ../../Infos.dat | awk '{ print $2 }'`
   actorb=`grep "Active_orbitals" ../../Infos.dat | awk '{ print $2 }'`
   inactorb=`grep "Inactive_orbitals" ../../Infos.dat | awk '{ print $2 }'`

   sed -i "s|nActEl=10 0 0|nActEl=$actelec 0 0|" ${newdir}.input
   sed -i "s|Ras2=10|Ras2=$actorb|" ${newdir}.input
   sed -i "s|Inactive=62|Inactive=$inactorb|" ${newdir}.input
fi

#
# Message for the user
#
cd ../../
cp $templatedir/ASEC/ASEC_direct_CASSCF.sh .
./update_infos.sh "Next_script" "ASEC_direct_CASSCF.sh" Infos.dat
echo ""
echo "***************************************************************************"
echo ""
echo "Run ASEC_direct_CASSCF.sh to generate the final coordinate file and submitt"
echo ""
echo "***************************************************************************"
echo ""

