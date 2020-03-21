#!/bin/bash
#
#
# Retrieving project name, parameter file, template folder and Tinker path from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
dowser=`grep "Dowser" Infos.dat | awk '{ print $2 }'`
#numatoms=`head -n2 Dynamic/$Project.gro | tail -n1 | awk '{ print $1 }'`

#
# Creating the directory where hydrogen atoms MM minimization will be done
# and putting the required files
#
if [[ -d Dynamic ]]; then
   echo ""
   echo "Dynamic folder already exists"
   echo " Aborting..."
   echo ""
   exit 0
fi

cp $templatedir/ASEC/update_infos.sh .

./update_infos.sh "******** ASEC " "********" Infos.dat 
./update_infos.sh "Step" 0 Infos.dat

mkdir Dynamic
cp $Project.pdb Dynamic
cd Dynamic
cp -r $templatedir/amber94.ff .
cd amber94.ff/
cp normalamino-h aminoacids.hdb
case $retstereo in
     AT) cp amino-rettrans aminoacids.rtp
     ;;
     11C) cp amino-ret11cis aminoacids.rtp
     ;;
     13C) cp amino-ret13cis aminoacids.rtp
     ;;
     nAT) cp amino-neutrans aminoacids.rtp
          cp neuamino-h aminoacids.hdb
     ;;
     new)
          for i in {1..54}; do
             lab=`head -n $i ../../new_charges | tail -n1 | awk '{ print $2 }'`
             charge=`head -n $i ../../new_charges | tail -n1 | awk '{ print $1 }'`
             sed -i "s/${lab}_RET/$charge/" amino-retnew
          done
          cp amino-retnew aminoacids.rtp
     ;;
esac
cd ..
cp $templatedir/residuetypes.dat .
cp $templatedir/standard-EM.mdp .
cp $templatedir/soglia .
#p $templatedir/pdb-to-gro.sh .

#cp $templatedir/carbret labelret

#
# pdb2gmx is the Gromacs utility for generating gro files and topologies
#
$gropath/gmx pdb2gmx -f $Project.pdb -o $Project.gro -p $Project.top -ff amber94 -water tip3p 2> grolog
checkgro=`grep 'Writing coordinate file...' grolog`
if [[ -z $checkgro ]]; then
   echo " An error occurred during the execution of pdb2gmx. Please look into grolog file"
   echo " No further operation performed. Aborting..."
   echo ""
   exit 0
else
   echo " new.gro and its topology were successfully generated"
   echo ""
   carica=`grep 'Total charge in system' grolog | awk '{ print $5 }'`
   carica=${carica/.*}
   if [[ $carica -ne 0 ]]; then
      if [[ $carica -gt 0 ]]; then
          poschg=0
          negchg=$carica
      else
          poschg=$(($carica*(-1)))
          negchg=0 
      fi
   else
         poschg=0
         negchg=0
   fi
   rm grolog
fi

echo ""
echo " $Project.gro was successfully generated from the optimized .pdb file"
echo ""

cd ..

#hess=b
#while [[ $hess != y && $hess != n ]]; do
#      echo ""
#      echo ""
#      echo " Do you want to update the Hessian along the iterative procedure? (y/n)"
#      echo ""
#      read hess
#done

hess=n

if [[ $hess ==  n ]]; then
   ./update_infos.sh "Update_hessian" "NO" Infos.dat
#   echo "Update_hessian NO" >> Infos.dat
fi

if [[ $hess ==  y ]]; then
   ./update_infos.sh "Update_hessian" "YES" Infos.dat
#   echo "Update_hessian YES" >> Infos.dat
   answer=0
   while  [[ $answer -ne 1 && $answer -ne 2 ]]; do
          echo ""
          echo " Select one option for updating the Hessian"
          echo ""
          echo " 1) From the previous step"
          echo " 2) From the QM/MM Optimization"
          echo ""
          read answer
   done
   if [ $answer -eq 1 ]; then
      ./update_infos.sh "Hessian_from" "previous" Infos.dat
#      echo "Hessian_from previous" >> Infos.dat
   else
      ./update_infos.sh "Hessian_from" "QMMM" Infos.dat
#      echo "Hessian_from QMMM" >> Infos.dat
   fi
fi

box=0
while  [[ $box -ne 1 && $box -ne 2 ]]; do
       echo ""
       echo " Now, select one of the following option"
       echo ""
       echo " 1) Continue the ASEC without solvent box"
       echo " 2) Embbed the protein in a Solvent box"
       echo ""
       read box
done

if [[ $box -eq 1 ]]; then
   cp $templatedir/ASEC/DynIt_list_chromo.sh .
   ./update_infos.sh "SolventBox" "NO" Infos.dat
#   echo "SolventBox NO" >> Infos.dat
   echo ""
   echo " Now run DynIt_list_chromo.sh to performe the MD"
   echo ""
else
   cp $templatedir/ASEC/Solvent_box.sh .
   ./update_infos.sh "SolventBox" "YES" Infos.dat
#   echo "SolventBox YES" >> Infos.dat
   echo ""
   echo " Now run Solvent_box.sh to generate the solvent box and thermalize it"
   echo ""
fi

