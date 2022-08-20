#!/bin/bash
#
# Retrieving the needed information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
uphess=`grep "Update_hessian" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
shell=`grep "Shell" Infos.dat | awk '{ print $2 }'`
chargechr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
fmnfad=`grep "Tail" Infos.dat | awk '{ print $2 }'`
redox=`grep "Redox" Infos.dat | awk '{ print $2 }'`

Molcami_OptSCF=`grep "Molcami_OptSCF" Infos.dat | awk '{ print $2 }'`

# Creating the directory calculations and copying all the needed files
# If calculations already exists, the script exits with a message
#
if [ -d calculations ]; then
   echo " Folder calculations already exists. Please check it and remove if necessary"
   echo " Terminating..."
   echo ""
   exit 0
fi

#
# Creating the tinker xyz and key files for the QM/MM geometry optimizations
#
echo "parameters $prm" > $Project.key
echo "parameters $prm" > ${Project}_nolink.key
echo "" >> $Project.key
echo "" >> ${Project}_nolink.key
chratoms=`head -n1 Chromophore/$chromophore.xyz | awk '{ print $1 }'`
#rm -f chromophore.xyz
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
if [[ $redox -eq 1 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/charges_FMN new_charges
fi
if [[ $redox -eq 2 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/charges_FMNH new_charges
fi
if [[ $redox -eq 3 && $fmnfad == "FMN" ]]; then
   cp $templatedir/ASEC/charges_FMNH2 new_charges
fi

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
# The charges of the MM atoms are equal to $carga (that is the meaning of the 0) to ensure 
# that the total charge of the QM+MM part is the desired charge.
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

mv ${Project}_nolink.key $Project-tk.key

#
# Adding the link atom to the tinker xyz and key files
#
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
# Reformating the tinket xyz file
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
# Collecting data to run the QM/MM geometry optimizations
#
mkdir calculations
newdir=${Project}_OptSCF
mkdir calculations/${newdir}

mv $Project.key calculations/${newdir}/${newdir}.key
cp $templatedir/ASEC/template_OptSCF calculations/${newdir}/template_OptSCF
if [[ $xx -gt 0 ]]; then
   sed -i "s/Charge = 0/Charge = $(($char+2))/g" calculations/${newdir}/template_OptSCF
else
   sed -i "s/Charge = 0/Charge = $(($chargechr+2))/g" calculations/${newdir}/template_OptSCF
fi
mv $Project-tk.xyz_2 calculations/${newdir}/${newdir}.xyz

#slurm
cp $templatedir/molcas.slurm.sh calculations/${newdir}/molcas-job.sh
#cp $templatedir/molcas-job.sh calculations/${newdir}/

cd calculations
cd ${newdir}/


NOME=${newdir}
NOME=${Project}_OptSCF
sed -i "s|NOMEPROGETTO|$NOME|" molcas-job.sh

no=$PWD
#sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|23000|" molcas-job.sh
sed -i "s|MEMORIA|20000|" molcas-job.sh
sed -i "s|hh:00:00|120:00:00|" molcas-job.sh

#
# Replacing PARAMETRI with current prm filename templateOPTSCF
#
sed -i "s|PARAMETRI|${prm}|" template_OptSCF
#sed -i "s/\*    rHidden = 4.0/    rHidden = 4.0/" template_OptSCF

sed -i "/export Project=/c\ export Project=${newdir}" molcas-job.sh

#
# Generation of the correct Molcas input from templateOPTSCF
#

if [[ $uphess == NO ]]; then
   sed -i "/RUNOLD/d" template_OptSCF
   sed -i "/OLDF/d" template_OptSCF
fi

mv template_OptSCF ${newdir}.input

#
# Message to the user
#
cd ../../
cp $templatedir/ASEC/ASEC.sh .
./update_infos.sh "Next_script" "ASEC.sh" Infos.dat
echo ""
echo "*************************************************************"
echo ""
echo "Run ASEC.sh to generate the final coordinate file and submitt"
echo ""
echo "*************************************************************"

