#!/bin/bash

#
# Retrieving information from Infos.dat
#
# If the ${Project}_VDZP_Opt exists, it means we are in calculations folder
#
if [ -d *_VDZP_Opt ]; then
   calculations="true"
   Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
   templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
   gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`
   prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
   tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
   dowser=`grep "Dowser" ../Infos.dat | awk '{ print $2 }'`
   chromophore=`grep "chromophore" ../Infos.dat | awk '{ print $2 }'`
   amber=`grep "AMBER" ../Infos.dat | awk '{ print $2 }'`
   chargechr=`grep "Chromo_Charge" ../Infos.dat | awk '{ print $2 }'`
   cola=`grep "Tail" ../Infos.dat | awk '{ print $2 }'`
else
   calculations="false"
   Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
   templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
   gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
   prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
   tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
   dowser=`grep "Dowser" Infos.dat | awk '{ print $2 }'`
   chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
   amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`
   chargechr=`grep "Chromo_Charge" Infos.dat | awk '{ print $2 }'`
   cola=`grep "Tail" Infos.dat | awk '{ print $2 }'`
fi
cd ESPF_charges

#
# Checking if the single point calculation ended succesfully
#
if [[ -f ${Project}_ESPF.out ]]; then
   if grep -q "Timing: Wall=" ${Project}_ESPF.out; then
      echo ""
      echo " Normal termination of ESPF charges."
      echo " Continuing ..."
      echo ""
   else
      echo ""
      echo " It seems to be that the ESPF calculation"
      echo " did not finish yet. Please check ..."
      echo ""
      exit 0
   fi
else
   echo ""
   echo " There is something wrong with the ESPF calculation."
   echo " Please check it. Finishing ..."
   echo ""
   exit 0
fi

#
# Colecting the charges from molcas output and rounding to the total charge
# to exactly the total charge of the QM part. 
# (there is alwais 0.002.. from molcas output)
#
numchr=`head -n1 ${chromophore}.xyz | awk '{ print $1 }'`
if [[ $cola == "FMN" ]]; then
   numqm=$(($numchr-20+1))
fi
if [[ $cola == "FAD" ]]; then
   numqm=$(($numchr-54+1))
fi

cat > new_charges.f << YOE
      Program new_gro
      implicit real*8 (a-h,o-z)
      character line80*80
      dimension charges($numqm),Ncharges($numqm),
     &          charges0($numqm)

      open(1,file='${Project}_ESPF.out',status='old')
      open(2,file='qm_charges',status='unknown')

CCC Reading the last obtained charges from molcas output

      do
         read(1,'(A80)',IOSTAT=io)line80
         if (io.gt.0) then
c            write(*,*) 'Check charges.out. Something was wrong'
            EXIT
         else if (io.eq.0) then
            k= index(line80,'Expectation values of the ESPF operators:')
            if (k.eq.7) then
               read(1,'(A80)')line80
               do i=1,$numqm
                  read(1,'(A35,f7.4)')line35,charges(i)
               enddo
            endif
         else if (io.lt.0) then
             goto 1
         endif
 1    enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c puting total charge exactly equal to 1

      itotal=$(($chargechr+2))*1000000
      do i=1,$numqm
         Ncharges(i)=nint(charges(i)*1000000)
      enddo
      k=sum(Ncharges)
c      write(*,*)k
 2    i=1
      do while (k.ne.itotal)
         if (k.gt.itotal) then
            Ncharges(i)=Ncharges(i)-1
         else
            Ncharges(i)=Ncharges(i)+1
         endif
         k=sum(Ncharges)
         i=i+1
      if (i.eq.$numqm+1) then
         goto 2
      endif
      enddo
      do i=1,$numqm
         charges0(i)=Ncharges(i)/1000000.0d0
         write(2,'(f9.6)')charges0(i)
      enddo
      end
YOE

gfortran new_charges.f -o new_charges.x
./new_charges.x

#
# Generating the whole charges of the chromophore + tail
# to be used by GROMACS
#

rm -f new_charges
touch new_charges
cont=1
for i in $(eval echo "{1..$numchr}"); do
   charge=`head -n $i template_charges | tail -n1 | awk '{ print $1 }'`
   atmtype=`head -n $(($i+2)) $chromophore.xyz | tail -n1 | awk '{ print $6 }'`
   if [[ $atmtype == "MM" ]]; then
      echo "$charge" >> new_charges
   fi
   if [[ $atmtype == "QM" ]]; then
      qmch=`head -n $cont qm_charges | tail -n1 | awk '{ print $1 }'`
      echo "$qmch" >> new_charges
      cont=$(($cont+1))
   fi
   if [[ $atmtype == "LQ" ]]; then
      qmch=`head -n $cont qm_charges | tail -n1 | awk '{ print $1 }'`
      echo "$qmch" >> new_charges
      cont=$(($cont+1))
   fi
   if [[ $atmtype == "LM" ]]; then
      qmch=`head -n $numqm qm_charges | tail -n1 | awk '{ print $1 }'`
      echo "$qmch" >> new_charges
   fi
#
#  1 for atoms of the tail as ASEC points
#  0 for the fixed atoms of the tail
#
   if [[ $atmtype == "XX" ]]; then
      echo "$charge" >> new_charges
   fi
   if [[ $atmtype == "FX" ]]; then
      echo "$charge" >> new_charges
   fi
done

#
# Creating the rtp file of the Chromophore + tail for GROMACS
#
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
      num=$numchr
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

#
# Message to user
#

if [[ $calculations == "false" ]]; then
   cp new_charges ../

   moldy="NVT"
#   while [[ $moldy != "NPT" && $moldy != "NVT" ]]; do
#      echo ""
#      echo " What ensemble will use for the MD along the iterative procedure (NPT or NVT)"
#      echo ""
#      echo ""
#      read moldy
#   done
   ../update_infos.sh "MD_ensemble" "$moldy" ../Infos.dat
   if [[ $moldy == "NVT" ]]; then
      echo ""
      echo "*************************************************************************"
      echo ""
      echo " Continuum with the MD_NPT.sh for equilibrating the volume of the system."
      echo ""
      echo "*************************************************************************"
      else
      echo ""
      echo " Continuum with the MD_NPT.sh."
      echo ""
      echo ""
   fi
   cp $templatedir/ASEC/MD_NPT.sh ../
   ../update_infos.sh "Next_script" "MD_NPT.sh" ../Infos.dat
else
   cp new_charges ../
   cp $templatedir/ASEC/Next_Iteration.sh ../../
   ../../update_infos.sh "Next_script" "Next_Iteration.sh" ../../Infos.dat
   echo ""
   echo "**********************************************************"
   echo " "
   echo " Go to the main folder and continue with: Next_Iteration.sh"
   echo ""
   echo "**********************************************************"
   echo ""
fi


