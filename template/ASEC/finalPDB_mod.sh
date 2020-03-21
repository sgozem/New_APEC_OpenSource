#!/bin/bash
#
# This script generates a PDB file of the final structure
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" ../Infos.dat | awk '{ print $2 }'`
solvent=`grep "SolventBox" ../Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" ../Infos.dat | awk '{ print $2 }'`
amber=`grep "AMBER" ../Infos.dat | awk '{ print $2 }'`

chr=`grep "chromophore" ../Infos.dat | awk '{ print $2 }'` 
numatoms=`grep "Shell" ../Infos.dat | awk '{ print $2 }'`
numchr=`head -n1 ../Chromophore/$chr.xyz | awk '{ print $1 }'`
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`

charge=`grep "Init_Charge" ../Infos.dat | awk '{ print $2 }'`
numions=${charge#-}  # module of charge

#
# Creating a directory where the PDB will be generated
# and copying the required files
#
echo " Preparing the final structure PDB..."
echo ""
folder=${Project}_VDZP_Opt

echo ""
echo " The current project is $Project. Checking the CAS/VDZP optimization..."
echo ""

#
# Checking if the optimization ended successfully, or if the assigned hours ended
#
if [ -d $folder ]; then
   cd $folder
   contr=`grep 'landing' $folder.out | awk '{ print $1 }'`
   if [[ $contr == "Happy" ]]; then
      echo " CAS/VDZP optimization ended successfully"
      echo ""
   else
      echo " The CAS/VDZP optimization did not finished well."
      echo " Please check if the assigned hours expired or if some error occurred."
      echo ""
      echo " If you need to restart, run restart.sh"
      echo " Otherwise find out what is the error and re-run sp_to_opt_VDZP_mod.sh."
      echo ""
      echo " 3rd_to_4th_mod.sh will terminate now!"
      echo ""
      exit 0
   fi
else
   echo " CAS/VDZP optimization folder was not found! Check what is wrong"
   echo " Terminating"
   exit 0
fi

#
# Retrieving the occupation numbers from CAS single point for checking purposes
#
grep -A2 "Natural orbitals and occupation numbers for root  1" ${folder}.out

#
# Asking the user if the occupation numbers are ok
#
echo ""
echo -n " Are all the occupation numbers in the range 0.02 to 1.98? (y/n)"
read answer
echo ""
contr=0
while [  $contr = 0 ]; do
      if [[ $answer != "y" && $answer != "n" ]]; then
         echo -n " Please answer y or n... -> "
         read answer
         echo ""
      else
         contr=1
         if [[ $answer == "y" ]]; then
            echo " Going ahead ..."
            echo ""
         else
            echo " You might have a problem with active space selection. To fix it:"
            echo " 1. - Use Molden to look at the orbitals in $folder.rasscf.molden"
            echo " 2. - Find the right orbital(s) to be placed in the active space"
            echo " 3. - Run the script alter_orbital_mod.sh"
            echo ""
            echo " Now 3rd_to_4th_mod.sh will terminate"
            echo ""
            cp $templatedir/ASEC/alter_orbital_mod.sh .
            exit 0
         fi
      fi
done
cd ..

mkdir ${Project}_finalPDB

#
# Getting the final xyz file from the $folder folder
# Terminates if the file is not found
#
if [ -f $folder/$folder.Final.xyz ]; then
   head -n $(($numatoms)) $folder/$folder.Final.xyz | tail -n $numchr > ${Project}_finalPDB/Chromo.xyz
   cp ../MD_ASEC/Best_Config_full.gro ${Project}_finalPDB
   cd ${Project}_finalPDB
   head -n1 Best_Config_full.gro > Protein.gro
   numfull=`head -n2 Best_Config_full.gro | tail -n1 | awk '{ print $1 }'`
   echo $(($numfull-$numchr-$numions)) >> Protein.gro
   grep "NA      NA\|CL      CL" Best_Config_full.gro > ions
   cp Best_Config_full.gro temp
   sed -i "/CHR /d" temp
   sed -i "/NA      NA/d" temp
   sed -i "/CL      CL/d" temp
   head -n $(($numfull-$numchr-$numions+3)) temp | tail -n $(($numfull-$numchr-$numions+1)) >> Protein.gro
   $gropath/gmx editconf -f Protein.gro -o Protein.pdb -label A
   sed -i "/ENDMDL/d" Protein.pdb
   cp ../../Chromophore/$chr.pdb .
cat > convert.f << YOE
      Program convert
      implicit real*8 (a-h,o-z)
      character label29*29,opls*8,line30*30

      open(1,file='$chr.pdb',status='old')
      open(2,file='Chromo.xyz',status='old')
      open(3,file='Chromo_new.pdb')
      do i=1,$numchr
         read(1,'(A)')label29
         read(2,'(12x,f11.6,1x,f11.6,1x,f11.6)')x,y,z
         write(3,'(A,f8.3,f8.3,f8.3)')label29,x,y,z
      enddo
      end
YOE
   gfortran convert.f -o convert.x
   ./convert.x
   cat Protein.pdb Chromo_new.pdb > ${Project}_new.pdb
   cp -r ../../Dynamic/$amber.ff .
   $gropath/gmx pdb2gmx -f ${Project}_new.pdb -o ${Project}_new.gro -p ${Project}_new.top -ff $amber -water tip3p 2> grolog

   if grep -q "NA      NA\|CL      CL" ions; then
      head -n1 ${Project}_new.gro > top
      echo "$numfull" >> top
      head -n $(($numfull-$numions+3)) ${Project}_new.gro | tail -n1 > vol
      head -n $(($numfull-$numions+2)) ${Project}_new.gro | tail -n $(($numfull-$numions)) > body
      cat top body ions vol > temp
   fi
   sed -i "s/HOH  /SOL  /g" temp
   mv temp ${Project}_new.gro

#
# Directions to update the parametrizations of the chromophore
#
   if [ -f ${Project}_new.gro ]; then
      cd ..
      cp $templatedir/ASEC/RESP_charges.sh .
     ../update_infos.sh "Next_script" "RESP_charges.sh" ../Infos.dat 
      echo ""
      echo "************************************************************************"
      echo ""
      echo " The files \"${Project}_new.gro\" and \"${Project}_new.pdb\" seems to be"
      echo " properly generated"
      echo ""
      echo " Continue with: RESP_charges.sh"
      echo ""
      echo "************************************************************************"
      echo ""  
   fi
else
   echo " $folder.xyz not found! Exiting... "
   echo ""
   exit 0
fi

