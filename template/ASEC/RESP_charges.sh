#!/bin/bash
#
# Reading data from Infos.dat
#

if [[ -d RESP_charges ]]; then
   echo " RESR_charges folder already exists "
   echo " Aborting ..."
   exit 0
fi

mkdir RESP_charges
cd RESP_charges

folder='VDZP_Opt'
Project=`grep "Project" ../../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../../Infos.dat | awk '{ print $2 }'`
chromo=`grep "CHR_RESP" ../../Infos.dat | awk '{ print $2 }'`
chr=`grep "Chromo_Charge" ../../Infos.dat | awk '{ print $2 }'`

cp ../${Project}_$folder/${Project}_$folder.Final.xyz final.xyz
cp ../${Project}_$folder/$prm.prm .
cp ../${Project}_$folder/${Project}_$folder.key key_file
cp $templatedir/ASEC/fitting_RESP.sh ../
cp ../../Dynamic/${Project}_box_sol.gro $Project.gro
cp ../../template_tk2gro .

#
# Here will be generated the gaussian input for computing the RESP charges (RET.com)
#


echo " $chr 1" > ${Project}_RESP.com
echo "New labels" > check_$chromo
echo "" >> check_$chromo

grep "QM \|MM " key_file | grep -v "QMMM" | awk '{ print $2 }' > qatoms
   
if [[ -f ../../chargefxx0 ]]; then
   awk '{ print $2 }' ../../chargefxx0 >> qatoms
fi

sort -n -o sorted qatoms
mv sorted qatoms
fin=`wc -l qatoms | awk '{ print $1 }'`
init=`head -n1 qatoms | awk '{ print $1 }'`
final=`head -n $fin qatoms | tail -n1 | awk '{ print $1 }'`

#  coordinates of the chromophore

count=0

for j in $(eval echo "{$init..$final}")
do
   count=$(($count+1))
   att=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $2 }' | awk '{print substr ($0, 0, 1)}'`
   xyz=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $3"   "$4"   "$5 }'`
   echo " $att    $xyz" >> ${Project}_RESP.com
   atnum=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $1 }'`
   atnumg=`head -n $(($j+1)) template_tk2gro | tail -n1 | awk '{ print $2 }'`
   labelg=`head -n $(($atnumg+2)) $Project.gro | tail -n1 | awk '{ print $2 }'`
   echo "$labelg" >> check_$chromo
done

echo "" >> ${Project}_RESP.com

diff -y $templatedir/ASEC/RESP/${chromo}_labels_order check_$chromo

echo ""
echo ""
echo " Please check that the labels of the Retinal (New labels)"
echo " are exactly in the same order as in the Reference."
echo " Otherwise, the retinal sequence atoms in the initial pdb"
echo " must be reorganized or the RESP input files (${chromo}-resp1.in"
echo " and ${chromo}-resp2.in) must be modified."
echo ""

order="b"
while [[ $order != "y" && $order != "n" ]]; do
   echo ""
   echo " Is the order the same? (y/n)"
   echo ""
   read order
done

if [[ $order == "n" ]]; then
   echo ""
   echo " Please modify the RESP input files ${chromo}-resp1.in and ${chromo}-resp2.in"
   echo " before continue"
   echo ""
   rm -r RESP_charges
   exit 0
else
   echo ""
   echo " Continuing ..."
   echo ""

fi

lines=`wc -l check_$chromo | awk '{ print $1 }'`
tail -n $(($lines-1)) check_$chromo > a
mv a check_$chromo

#
#taking the charges from parameters file
#

ncharges=`grep -c "charge " $prm.prm`
grep "charge " $prm.prm > charges

# writing the point charges to RET.com, not including the retinal and puting to ZERO the charge of CD (close to LAH): Init[1]-2

total=`head -n1 final.xyz | awk '{ print $1 }'`

rm -f tempiLOV
echo "" >> final.xyz

cat > charges.f << YOE
      Program charges
      implicit real*8 (a-h,o-z)
      character line3*3,line7*7
      dimension coorx($total),coory($total),coorz($total)
     &          ,charge(9999),indi($total)

      open(1,file='final.xyz',status='old')
      open(2,file='charges',status='old')
      open(3,file='tempCHR',status='unknown')

      read(1,*)
      do i=1,$total
         read(1,*)ini,line3,coorx(i),coory(i),coorz(i),indi(i)
      enddo

      do i=1,$ncharges
         read(2,*)line7,ind,value
         charge(ind)=value
      enddo

      do i=1,$total
         write(3,'(1x,3(f10.6,3x),f9.6)')coorx(i),coory(i),
     &   coorz(i),charge(indi(i))
      enddo
      end
YOE

gfortran charges.f -o charges.x
./charges.x
rm charges.x

sed -e "${init},${final}s/.*/DELETE/" tempCHR > a
mv a tempCHR

sed -i "/DELETE/d" tempCHR

cat ${Project}_RESP.com tempCHR > a
mv a ${Project}_RESP.com

echo "" >> ${Project}_RESP.com 

rm charges

state="22"
while [[ $state -ne 0 && $state -ne 1 ]]; do
   echo ""
   echo " Please select what electronic state do you want to parametrize?"
   echo ""
   echo "0: Ground state"
   echo "1: First excited state"
   read state
done

if [[ $state -eq 0 ]]; then
   cp $templatedir/ASEC/RESP/RESP.com .
   ../../update_infos.sh "State" 0 ../../Infos.dat
else
   cp $templatedir/ASEC/RESP/EXCITED.com RESP.com
   ../../update_infos.sh "State" 1 ../../Infos.dat
fi

cat RESP.com ${Project}_RESP.com > temp10 
mv temp10 ${Project}_RESP.com
rm RESP.com

cp $templatedir/ASEC/gaussian.sh .
sed -i "s/PROJECTO/${Project}_RESP/" gaussian.sh

#slurm
#sbatch gaussian.sh
qsub gaussian.sh

cd ..
cp $templatedir/ASEC/fitting_RESP.sh .

sed -i "s/STATFIT/TNK/" fitting_RESP.sh
echo ""
echo ""
echo " Gaussian calculation submitted (HF/6-31G* plus ASEC) for generating the"
echo " Molecular Electrostatic Potential (MEP)."
echo ""
echo " When this calculation is done, execute fitting_RESP.sh"
echo ""
echo ""

