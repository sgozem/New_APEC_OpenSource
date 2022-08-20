#!/bin/bash
#
# Reading project name, parm file and templatedir from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
multichain=`grep "MultChain" Infos.dat | awk '{ print $2 }'`
step=`grep "Step" Infos.dat | awk '{ print $2 }'`
#solventbox=`grep "SolventBox" Infos.dat | awk '{print $2}'`

#cp $Project.ndx Dynamic

cd Dynamic/

if [[ $step -eq 0 ]]; then
#   if [[ $solventbox == YES ]]; then
#      mv $Project.gro ${Project}_no_Sol.gro
#      mv $Project.top ${Project}_no_Sol.top
#      cp SolventBox/output/final-${Project}_box_sol.gro $Project.gro
#      cp SolventBox/$Project.top . 
#   fi
#
# Generation of the .ndx file. This wil be used along the iterative procedure
#

   ../update_infos.sh "MDRelax" "Cavity" ../Infos.dat
   cp $templatedir/ASEC/ndx-maker_mod.sh .

   backb=5
   while  [[ $backb -ne 0 && $backb -ne 1 ]]; do
          echo ""
          echo " Please type 1 if you want to relax the backbone, 0 otherwise"
          echo ""
          echo ""
          read backb
   done

   answer=0
   while  [[ $answer -ne 1 && $answer -ne 2 && $answer -ne 3 && $answer -ne 4 ]]; do
          echo ""
          echo " Now, select one option for determining the side chains of the cavity:"
          echo ""
          echo " 1) Use a previously defined cavity files"
          echo " 2) Generate the .ndx for the 4 A cavity"
          echo " 3) Generate the .ndx for the 6 A cavity"
          echo " 4) Generate the .ndx for the 8 A cavity"
          echo ""
          read answer
   done

   case $answer in
        1) 
#           backb=0
           ../update_infos.sh "BackBoneMD" $backb ../Infos.dat
           ../update_infos.sh "RadiusMD" "ND" ../Infos.dat
           echo  " Enter the full path, including the name, of the cavity file"
           read CavityFile
           cp $CavityFile ../cavity
           ../update_infos.sh "CavityFile" "YES" ../Infos.dat

           waters=10
           while  [[ $waters -ne 4 && $waters -ne 6 && $waters -ne 8 ]]; do
              echo  ""
              echo  ""
              echo " Sometimes the water molecules are not included in the cavity file."
              echo " Anyway, select the size of the cavity for selecting the water molecules (4, 6 or 8 A)" 
              echo ""
           echo ""
           read waters
           done

           ../update_infos.sh "Waters" $waters ../Infos.dat
#           echo "Waters $waters" >> ../Infos.dat

           echo "" 
           echo '2' > choices.txt
           echo 'q' >> choices.txt
           $gropath/gmx make_ndx -f $Project.gro -o $Project.ndx < choices.txt
           rm choices.txt
           ./ndx-maker_mod.sh $Project 2 $backb
        ;;
        2) 
#           backb=0
           ../update_infos.sh "BackBoneMD" $backb ../Infos.dat
           ../update_infos.sh "RadiusMD" 4 ../Infos.dat
           echo "" 
           echo '2' > choices.txt
           echo 'q' >> choices.txt
           $gropath/gmx make_ndx -f $Project.gro -o $Project.ndx < choices.txt
           rm choices.txt
           ./ndx-maker_mod.sh $Project 2 $backb
         ;;
        3) 
#           backb=0
           ../update_infos.sh "BackBoneMD" $backb ../Infos.dat
           ../update_infos.sh "RadiusMD" 6 ../Infos.dat
           echo "" 
           echo '2' > choices.txt
           echo 'q' >> choices.txt
           $gropath/gmx make_ndx -f $Project.gro -o $Project.ndx < choices.txt
           rm choices.txt
           ./ndx-maker_mod.sh $Project 2 $backb
        ;;
        4)
#           backb=0
           ../update_infos.sh "BackBoneMD" $backb ../Infos.dat
           ../update_infos.sh "RadiusMD" 8 ../Infos.dat
           echo "" 
           echo '2' > choices.txt
           echo 'q' >> choices.txt
           $gropath/gmx make_ndx -f $Project.gro -o $Project.ndx < choices.txt
           rm choices.txt
           ./ndx-maker_mod.sh $Project 2 $backb
        ;;
   esac
   cp $Project.ndx ../
fi

#cp $templatedir/update_infos.sh ../

cp $templatedir/ASEC/moldynamic_parall.sh .
sed -i "s/update_infos.sh \$checkrest/update_infos.sh/g" moldynamic_parall.sh

./moldynamic_parall.sh

cd ..
cp $templatedir/ASEC/MD_ASEC.sh .
cp $templatedir/Analysis_MD.sh .

echo ""
echo " Wait for molecular dynamics to end, then run MD_ASEC.sh "
echo ""

