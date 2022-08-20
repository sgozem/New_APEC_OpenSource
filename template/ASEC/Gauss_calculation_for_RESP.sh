#!/bin/bash
#
# Reading data from Infos.dat
#

echo ""
echo ""
echo " Enter the name of the Gaussian calculation folder"
echo ""
read gaussname

if [[ -d $gaussname ]]; then
   echo " $gausname folder already exists "
   echo " Aborting ..."
   exit 0
fi

echo " please wait ..."

mkdir $gaussname
cd $gaussname

folder="VDZP_Opt"
Project=`grep "Project" ../../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../../Infos.dat | awk '{ print $2 }'`
iLOV="iLOV_chain"
chr=`grep "Chromo_Charge" ../../Infos.dat | awk '{ print $2 }'`

cp ../${Project}_$folder/${Project}_$folder.xyz final.xyz
cp ../${Project}_$folder/$prm.prm .
cp ../${Project}_$folder/${Project}_$folder.key key_file
cp ../../Dynamic/${Project}_box_sol.gro $Project.gro
cp ../../template_tk2gro .

#
# Here will be generated the gaussian input for computing gaussian calculations
#

echo " $chr 1" > ${Project}_GAUSS.com
echo ""

if [[ $iLOV == "iLOV" ]]; then

   interv_qm=`grep -w -c "QM" key_file`

   for i in $(eval echo "{1..$interv_qm}")
   do
      init[$i]=`grep -w -m $i "QM" key_file | tail -n1 | awk '{print ((-1*$2))}'`
      final[$i]=`grep -w -m $i "QM" key_file | tail -n1 | awk '{print $3}'`
   done
fi
if [[ $iLOV == "iLOV_chain" ]]; then
   grep "QM \|MM " key_file | grep -v "QMMM" | awk '{ print $2 }' > qatoms
   sort -n -o sorted qatoms
   mv sorted qatoms
   fin=`wc -l qatoms | awk '{ print $1 }'`
   interv_qm=1
   init[1]=`head -n1 qatoms | awk '{ print $1 }'`
   final[1]=`head -n $fin qatoms | tail -n1 | awk '{ print $1 }'`
fi

   #  coordinates of the chromophore

count=0

for i in $(eval echo "{1..$interv_qm}")
do
   for j in $(eval echo "{${init[$i]}..${final[$i]}}")
   do
      count=$(($count+1))

      att=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $2 }' | awk '{print substr ($0, 0, 1)}'`
      xyz=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $3"   "$4"   "$5 }'`
      echo " $att    $xyz" >> ${Project}_GAUSS.com
      atnum=`head -n $(($j+1)) final.xyz | tail -n1 | awk '{ print $1 }'`
      atnumg=`head -n $(($j+1)) template_tk2gro | tail -n1 | awk '{ print $2 }'`
      labelg=`head -n $(($atnumg+2)) $Project.gro | tail -n1 | awk '{ print $2 }'`
   done
done

echo "" >> ${Project}_GAUSS.com


#taking the charges from parameters file
ncharges=`grep -c "charge " $prm.prm`
grep "charge " $prm.prm > charges

# writing the point charges to RET.com, not including the retinal and puting to ZERO the charge of CD (close to LAH): Init[1]-2

total=`head -n1 final.xyz | awk '{ print $1 }'`

if [ -f tempiLOV ]; then
   rm tempiLOV
fi
echo "" >> final.xyz

#CD=`grep " LAH \| HLA " final.xyz | awk '{ print $7 }'`
#LAH=`grep " LAH \| HLA " final.xyz | awk '{ print $1 }'`

cat > charges.f << YOE
      Program charges
      implicit real*8 (a-h,o-z)
      character line3*3,line7*7
      dimension coorx($total),coory($total),coorz($total)
     &          ,charge(9999),indi($total)

      open(1,file='final.xyz',status='old')
      open(2,file='charges',status='old')
      open(3,file='tempiLOV',status='unknown')

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

for i in $(eval echo "{1..$interv_qm}"); do
   sed -e "${init[$i]},${final[$i]}s/.*/DELETE/" tempiLOV > a
   mv a tempiLOV
done

sed -i "/DELETE/d" tempiLOV

cat ${Project}_GAUSS.com tempiLOV > a
mv a ${Project}_GAUSS.com

echo "" >> ${Project}_GAUSS.com 

rm charges

cp $templatedir/ASEC/GAUSS_CALC.com .

cat GAUSS_CALC.com ${Project}_GAUSS.com > temp10 
mv temp10 ${Project}_GAUSS.com
rm GAUSS_CALC.com

cp $templatedir/ASEC/gaussian.sh .
sed -i "s/PROJECTO/${Project}_GAUSS/" gaussian.sh
sed -i "s/SBATCH -t 9:00:00/SBATCH -t 220:00:00/" gaussian.sh
sed -i "s/SBATCH --mem=12000MB/SBATCH --mem=32000MB/" gaussian.sh

#slurm
#sbatch gaussian.sh
#qsub gaussian.sh

cd ..

echo ""
echo ""
echo " Gaussian input calculation created into $gaussname ..." 
echo " Just define the level of calculation and basis set and submit the job"
echo ""
echo ""
echo ""

