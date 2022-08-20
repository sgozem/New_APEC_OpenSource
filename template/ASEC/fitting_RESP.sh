#!/bin/bash
#
# Reading data from Infos.dat
#

if [[ -f Infos.dat ]]; then
   RESP="GRO"
else
   RESP="TNK"
fi

if [[ $RESP != "GRO" && $RESP != "TNK" ]]; then
   echo " Please define RESP type in this script."
   echo " Aborting ..."
   exit 0
fi

if [[ $RESP == "GRO" ]]; then
   Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
   templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
   chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
   Chromo=`grep "CHR_RESP" Infos.dat | awk '{ print $2 }'`
   dimer=`grep "Dimer" Infos.dat | awk '{ print $2 }'`
   chr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
else
   Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
   templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
   chromophore=`grep "chromophore" ../Infos.dat | awk '{ print $2 }'`
   Chromo=`grep "CHR_RESP" ../Infos.dat | awk '{ print $2 }'`
   cp ../Chromophore/$chromophore.xyz RESP_charges
   dimer=`grep "Dimer" ../Infos.dat | awk '{ print $2 }'`
   chr=`grep "Chromo_Charge" ../Infos.dat | awk '{ print $2 }'`
fi

cd RESP_charges

if [[ $dimer == "NO" ]]; then
   if [[ -f ${Project}_RESP.log ]]; then
      cp ${Project}_RESP.log ${Project}_RESP.out
      if grep -q "Normal termination of Gaussian" ${Project}_RESP.out; then
         echo ""
         echo " Normal termination of Gaussian."
         echo " Continuing ..."
         echo ""
      else
         echo ""
         echo " There is something wrong with the Gaussian calculation."
         echo " Please check it. Finishing ..."
         echo ""
         exit 0
      fi
   else
      echo ""
      echo " It seems to be that the Gaussian calculation"
      echo " did not finish yet. Please check ..."
      echo ""
      exit 0
   fi

   cp $templatedir/ASEC/RESP/readit.f .
   cp $templatedir/ASEC/RESP/$Chromo-resp1.in .
   cp $templatedir/ASEC/RESP/$Chromo-resp2.in .

   sed -i "s/CHARGE/$chr/g" $Chromo-resp1.in
   sed -i "s/CHARGE/$chr/g" $Chromo-resp2.in

   grep "Atomic Center " ${Project}_RESP.out > a
   grep "ESP Fit" ${Project}_RESP.out > b
   grep "Fit    " ${Project}_RESP.out > c

   i=`grep -c "Atomic Center " ${Project}_RESP.out`
   j=`grep -c "ESP Fit" ${Project}_RESP.out`
   sed -i "s/num1/$i/" readit.f
   sed -i "s/num2/$j/" readit.f

   if [ -f esp.dat ]; then
      rm esp.dat
   fi

   gfortran readit.f
   ./a.out
   rm -f a b c a.out readit.o

   $templatedir/ASEC/RESP/resp -O -i $Chromo-resp1.in -o $Chromo-resp1.out -p $Chromo-resp1.pch -t $Chromo-resp1.chg -e esp.dat
   $templatedir/ASEC/RESP/resp -O -i $Chromo-resp2.in -o $Chromo-resp2.out -p $Chromo-resp2.pch -t $Chromo-resp2.chg -e esp.dat -q $Chromo-resp1.chg

   lines=`wc -l $Chromo-resp2.chg | awk '{ print $1 }'` 

   if [ -f new_charges ]; then
      rm new_charges
   fi

   for i in $(eval echo "{1..$lines}"); do
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $1 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $2 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $3 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $4 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $5 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $6 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $7 }'`
      echo "$char" >> new_charges
      char=`head -n $i $Chromo-resp2.chg | tail -n1 | awk '{ print $8 }'`
      echo "$char" >> new_charges
   done

   #paste new_charges check_RET | awk '{ print $1"   "$2 }' > a
   #mv a new_charges
   cp new_charges ../
fi

numatm=`head -n1 $chromophore.xyz | awk '{ print $1 }'`
end=`grep -n "End\|end\|END" $chromophore.xyz | cut -f1 -d:`
cat > write_charges.f << YOE
      Program write_charges
      implicit real*8 (a-h,o-z)
      character label*3,opls*8,line30*30
c      dimension charge($numatm)

      open(1,file='$chromophore.xyz',status='old')
      open(2,file='new_charges',status='old')
      open(3,file='new_rtp',status='unknown')

CCCCCCCCC Number of atoms of the solute
      num=$numatm
CCCCCCCCC
      read(1,*)
      read(1,*)
      write(3,'(A)')'[ CHR ]'
      write(3,'(A)')' [ atoms ]'
      do i=1,num
         read(1,*)label,x,y,z,opls
         read(2,*)charge
         write(3,'(1x,A,4x,A,7x,f10.6,3x,i5)')label,opls,charge,i
      enddo
      write(3,'("[bonds]")')
      do i=num+3,$end-1
         read(1,'(A)')line30
         write(3,'(A)')line30
      enddo
      write(3,*)
      end
YOE

gfortran write_charges.f -o write_charges.x
./write_charges.x

if [[ $dimer == "NO" ]]; then
   if grep -q "Statistics of the fitting" $Chromo-resp2.out; then
      echo ""
      echo ""
      echo " RESP charges seems to be properly computed!"
      echo ""
      if [[ $RESP == "TNK" ]]; then
         echo " Go to the main folder (cd ..) and continue with \"Next_Iteration.sh\""
         echo ""
         cp $templatedir/ASEC/Next_Iteration.sh ../../
      else
         echo " Continue with \"Solvent_box.sh\""
         cp $templatedir/ASEC/Solvent_box.sh ../
         echo ""
      fi
   else
      echo ""
      echo ""
      echo " There is something wrong with the RESP fitting"
      echo ""
   fi
else
   if [[ $RESP == "TNK" ]]; then
      echo " Go to the main folder (cd ..) and continue with \"Next_Iteration.sh\""
      echo ""
      cp $templatedir/ASEC/Next_Iteration.sh ../../
   else
      echo " Continue with \"Solvent_box.sh\""
      cp $templatedir/ASEC/Solvent_box.sh ../
      echo ""
   fi
fi


