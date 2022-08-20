#!/bin/bash
#
# Reading information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
step=`grep "Step" Infos.dat | awk '{ print $2 }'`
heat=`grep "HeatMD" Infos.dat | awk '{ print $2 }'`
equi=`grep "EquiMD" Infos.dat | awk '{ print $2 }'`
prod=`grep "ProdMD" Infos.dat | awk '{ print $2 }'`
#charge=`grep "Init_Charge" Infos.dat | awk '{ print $2 }'`
numparallel=`grep "Parallel_MD" Infos.dat | awk '{ print $2 }'`
dimer=`grep "Dimer" Infos.dat | awk '{ print $2 }'`
chromophore=`grep "chromophore" Infos.dat | awk '{ print $2 }'`
nna=`grep "Added_NAs" Infos.dat | awk '{ print $2 }'`
ncl=`grep "Added_CLs" Infos.dat | awk '{ print $2 }'`

confs=100
if [[ $numparallel -gt 1 ]]; then
   for i in $(eval echo "{1..$numparallel}")
   do
      if [[ -d Dynamic/seed_$i/output ]]; then
         echo ""
         echo " *********************************************************************"
         echo "                      Warning!"
         echo ""
         echo " MD seed_$i/output directory already excist. We are goint to use it..."
         echo ""
         echo " *********************************************************************"
         echo ""
      else
         echo ""
         echo " Collecting data from the iRODS, please wait ..."
         echo ""
         seeddir=`grep "MD_seed_$i" Infos.dat | awk '{ print $2 }'`
         dir=`basename $seeddir`
         iget -r /arctic/projects/CHEM9C4/$USER/$dir Dynamic/seed_$i
         if [[ -f Dynamic/seed_$i/$dir/md.log ]]; then
            mv Dynamic/seed_$i/$dir Dynamic/seed_$i/output
            irm -r /arctic/projects/CHEM9C4/$USER/$dir
         else
            echo ""
            echo "************************************************************************"
            echo ""
            echo " It seems that the MD is still running or it did not finish properly"
            echo ""
            echo "************************************************************************"
            echo ""
            exit 0
         fi
      fi
   done
else
   if [[ -d Dynamic/output ]]; then
      echo ""
      echo " *********************************************************************"
      echo "                      Warning!"
      echo ""
      echo " MD output directory already excist. We are goint to use it..."
      echo ""
      echo " *********************************************************************"
      echo ""
   else
      echo ""
      echo " Collecting data from the iRODS, please wait ..."
      echo ""
      tempdir=`grep "tempdir" Infos.dat | awk '{ print $2 }'`
      dir=`basename $tempdir`
      iget -r /arctic/projects/CHEM9C4/$USER/$dir Dynamic
      if [[ -f Dynamic/$dir/md.log ]]; then
         mv Dynamic/$dir Dynamic/output
         irm -r /arctic/projects/CHEM9C4/$USER/$dir
      else
         echo ""
         echo "************************************************************************"
         echo ""
         echo " It seems that the MD is still running or it did not finish properly"
         echo ""
         echo "************************************************************************"
         echo ""
         exit 0
      fi
   fi
fi

#
# Reading fom the MD output (tpr) to select 100 equally separated configurations to
# generate the ASEC configuration.
#
# -"paso" indicates the time step in ps to write configurations in the ASEC configuration
# -"skip" is the number configurations between the selected ones. It is based
#  on the number of confugurations written in the tpr file
# -"init" is the time in ps of the fist selected configuration
#
wrt=`grep "nstxtcout" Dynamic/dynamic_sol_NVT.mdp | awk '{ print $3 }'`
dt=`grep "dt                      =" Dynamic/dynamic_sol_NVT.mdp | awk '{ print $3 }'`
paso=$(echo "scale=0; ($dt*$wrt)/1" | bc)

mkdir MD_ASEC
cd MD_ASEC

# w means written configurations
heatw=$(echo "scale=0; $heat/$paso" | bc)
equiw=$(echo "scale=0; $equi/$paso" | bc)
prodw=$(echo "scale=0; $prod/$paso" | bc)

echo " These are the number of writen configurations:"
echo "heat  equi  prod"
echo "$heatw $equiw $prodw"
echo " Writing step $paso ps"

#if [[ $heat -eq 0 ]]; then
#   skip=$((($prodw-(300/$paso))/$confs))
#   init=$(($prod-300+$skip*$paso))
#else
   skip=$(($prodw/$confs))
   init=$(($skip*$paso+$heat+$equi))
#fi
echo "init $init" ps
echo "skip $skip" configurations

#
# After selecting the configurations using "trjconv", "trjorder" is used 
# to order the water molecules by distance from the chromophore. 
# This is needed in the future to select the 22 A shell to run the QMMM calculations
#
if [[ $numparallel -gt 1 ]]; then
   for i in $(eval echo "{1..$numparallel}")
   do
      $gropath/gmx trjconv -s ../Dynamic/seed_$i/output/${Project}_box_sol.tpr -f ../Dynamic/seed_$i/output/${Project}_box_sol.xtc -b $init -skip $skip -o Selected_100_seed_$i.gro << EOF
0
EOF
      $gropath/gmx trjorder -f Selected_100_seed_$i.gro -s ../Dynamic/seed_$i/output/${Project}_box_sol.tpr -n ../Dynamic/seed_$i/output/${Project}_box_sol.ndx -o Ordered_100_seed_$i.gro
# << EOF
#16
#13
#EOF
      if [[ $i -eq 1 ]]; then
         cp Selected_100_seed_$i.gro Selected_100.gro
         cp Ordered_100_seed_$i.gro Ordered_100.gro
      else
         cat Selected_100.gro Selected_100_seed_$i.gro > temp
         cat Ordered_100.gro Ordered_100_seed_$i.gro > temp2 
         mv temp Selected_100.gro
         mv temp2 Ordered_100.gro
      fi
   done
else
   $gropath/gmx trjconv -s ../Dynamic/output/${Project}_box_sol.tpr -f ../Dynamic/output/${Project}_box_sol.xtc -b $init -skip $skip -o Selected_100.gro << EOF
0
EOF
   $gropath/gmx trjorder -f Selected_100.gro -s ../Dynamic/output/${Project}_box_sol.tpr -n ../Dynamic/output/${Project}_box_sol.ndx -o Ordered_100.gro
# << EOF
#16
#13
#EOF
fi

#if [[ $charge -ne 0 ]]; then
#   if [[ $charge -lt 0 ]]; then
#      nions=$(echo "-1*$charge" | bc)
#   else
#      nions=$charge
#   fi
#else
#   nions=0
#fi

nions=$(($nna+$ncl))

echo ""
echo " wait ..."
echo ""

if [[ $step -eq 0 ]]; then

#
# The first configuration out of the 100 selected will be used to select the amount of atoms contained 
# in the 22 A shell
#

   att=`head -n2 Ordered_100.gro | tail -n1 | awk '{ print $1 }'`
   head -n $(($att+3)) Ordered_100.gro > first.gro
   att3=$(($att+3))
   att1=$(($att+1))
   vol=`tail -n1 first.gro | awk '{ print $0 }'`

   echo -e "mol new first.gro type gro" > selection.tcl
   echo -e "mol delrep 0 top" >> selection.tcl
#
# just selecting the environment of the first CHR of the dimer
#
   if [[ $dimer == "YES" ]]; then
      grep -n 'CHR ' first.gro | cut -d : -f 1 > temporal
      awk '$1=$1-2' temporal > t
      conta=`grep -c "CHR " first.gro | awk '{ print $1 }'`
      head -n $(($conta/2)) t > temporal1
      tail -n $(($conta/2)) t > temporal2
      sel1=`cat temporal1 | tr "\n" " " | awk '{ print $0 }'`
      sel2=`cat temporal2 | tr "\n" " " | awk '{ print $0 }'`
#
#  There are cases where some residues are not within the 22 A. So, protein has been added. In addition, for some reason, the terminal residue is not selected when "protein" is used,
#  so "same residue as name OC1 ..." was used to ensure that the whole protein is selected  
#
      var1="set selec [ atomselect top \"protein or (same residue as name OC1 and same residue as name OC2) or same residue as (all within 22 of (serial $sel1)) not (resname NA CL) and not (serial $sel2) \" ]"      
   fi
   if [[ $dimer == "NO" ]]; then
#
#  The selection of the shell is just based on the QM + MM atoms, because sometimes the tail can be to long
#
       rm -f sele
       grep -n 'CHR ' first.gro | cut -d : -f 1 > temporal
       awk '$1=$1-2' temporal > t
       mv t temporal
       chratoms=`head -n1 ../Chromophore/$chromophore.xyz | awk '{ print $1 }'`
       for i in $(eval echo "{1..$chratoms}"); do
          atmtype=`head -n $(($i+2)) ../Chromophore/$chromophore.xyz | tail -n1 | awk '{ print $6 }'`
          if [[ $atmtype == "MM" || $atmtype == "QM" || $atmtype == "LM" || $atmtype == "LQ" ]]; then
             head -n $i temporal | tail -n1 | awk '{ print $1 }' >> sele
          fi
       done
       sel=`cat sele | tr "\n" " " | awk '{ print $0 }'`
       var1="set selec [ atomselect top \"protein or (same residue as name OC1 and same residue as name OC2) or same residue as (all within 22 of (serial $sel)) and not (resname NA CL) \" ]"
   fi
   echo -e "$var1" >> selection.tcl
   echo -e 'set fileId [open tempnum "w"]' >> selection.tcl
   echo -e 'set numat [$selec num]' >> selection.tcl
   echo -e 'puts $fileId $numat' >> selection.tcl
   echo -e 'close $fileId' >> selection.tcl
   echo -e "exit" >> selection.tcl
   vmd -e selection.tcl -dispdev text
   shel=`head -n1 tempnum | awk '{ print $1 }'`
   if [[ $dimer == "YES" ]]; then
      shell=$(($shel+$nions+$conta/2))
   fi
   if [[ $dimer == "NO" ]]; then
      shell=$(($shel+$nions))
   fi
   ../update_infos.sh "Shell" "$shell" ../Infos.dat
   rm first.gro selection.tcl tempnum

#
# The Templates gro2tk and tk2gro will be modified to add the ions, the CHR and the waters
#
   cd ..
   cp Templates/template_* .
   numchrom=`grep -c "CHR " Minimize_$Project/$Project.gro`
   protatoms=`head -n2 Minimize_$Project/$Project.gro | tail -n1 | awk '{ print $1 }'`
   wat=`grep "DOWSER_wat" Infos.dat | awk '{ print $2 }'`
   nwat=$(($wat+$wat+$wat))

   if [[ $protatoms -gt $shell ]]; then
      echo "##################################################################"
      echo ""
      echo " The number of atoms of the protein is larger that the number"
      echo " of atoms of the shell."
      echo " Please, increase the distance criterion defining the shell"
      echo " but, the total number of atoms in the final ASEC configuration"
      echo " must be lower than 1000000."
      echo ""
      echo " Aborting ..."
      echo ""
      echo "#################################################################"
      exit 0
   fi
#
# when solvent box is used the number of atoms of the protein (protatoms) is only used by 
# MD_ASEC.f to computing the rms
#
   ./update_infos.sh "protatoms" $protatoms Infos.dat

   echo ""
   echo " wait ..."
   echo ""

#
# The waters molecules kept by dowser were considered in the templates before CHR
# because there was no CHR at that time. So, they will be removed now.
#
   list=`head -n2 Templates/$Project.gro | tail -n1 | awk '{ print $1 }'`
   head -n $(($list+1-$nwat)) template_gro2tk > temp1
   mv temp1 template_gro2tk
   head -n $(($list+1-$nwat)) template_tk2gro > temp1
   mv temp1 template_tk2gro
   
#
# The template files, which convert gro to tk and viceversa will be modified in order to add
# the water molecules and the ions 
#
#CHROMOPHORE block
   cont=1
   for i in $(eval echo "{$(($protatoms-$nwat-$numchrom+1))..$((protatoms-$nwat))}"); do
      if [[ $i -le 9999 ]]; then
         echo " $i    $(($shell-$numchrom+$cont))" >> template_gro2tk
      else
         echo "$i    $(($shell-$numchrom+$cont))" >> template_gro2tk
      fi
      cont=$(($cont+1))
   done
#WATER block
   cont=1
   for i in $(eval echo "{$(($protatoms-$nwat+1))..$(($shell-$nions))}"); do
      if [[ $i -le 9999 ]]; then
         echo " $i    $(($protatoms-$nwat-$numchrom+$cont))" >> template_gro2tk
      else
         echo "$i    $(($protatoms-$nwat-$numchrom+$cont))" >> template_gro2tk
      fi
      cont=$(($cont+1))
   done
#IONS block
   if [[ $nions -gt 0 ]]; then
      cont=1
      for i in $(eval echo "{$(($shell-$nions+1))..$shell}"); do
         if [[ $i -le 9999 ]]; then
            echo " $i    $(($shell-$numchrom-$nions+$cont))" >> template_gro2tk
         else
            echo "$i    $(($shell-$numchrom-$nions+$cont))" >> template_gro2tk
         fi
         cont=$(($cont+1))
      done
   fi

#################################################################################################

#WATER block
   cont=1
   for i in $(eval echo "{$(($protatoms-$nwat-$numchrom+1))..$(($shell-$numchrom-$nions))}"); do
      if [[ $i -le 9999 ]]; then
         echo " $i    $(($protatoms-$nwat+$cont))" >> template_tk2gro
      else
         echo "$i    $(($protatoms-$nwat+$cont))" >> template_tk2gro
      fi
      cont=$(($cont+1))
   done

#IONS block
   if [[ $nions -gt 0 ]]; then
      cont=1
      for i in $(eval echo "{$(($shell-$numchrom-$nions+1))..$(($shell-$numchrom))}"); do
         if [[ $i -le 9999 ]]; then
            echo " $i    $(($shell-$nions+$cont))" >> template_tk2gro
         else
            echo "$i    $(($shell-$nions+$cont))" >> template_tk2gro
         fi
         cont=$(($cont+1))
      done
   fi
#CHROMOPHORE block
   cont=1
   for i in $(eval echo "{$(($shell-$numchrom+1))..$shell}"); do
      if [[ $i -le 9999 ]]; then
         echo " $i    $(($protatoms-$nwat-$numchrom+$cont))" >> template_tk2gro
      else
         echo "$i    $(($protatoms-$nwat-$numchrom+$cont))" >> template_tk2gro
      fi
      cont=$(($cont+1))
   done

else
#
# Due to the link atom added in Step_0, shell was increased
#
   shell=`grep "Shell" ../Infos.dat | awk '{ print $2 }'`
#   shell=$(($shell-1))
#   ../update_infos.sh "Shell" "$shell" ../Infos.dat
   cd ..
   wat=`grep "DOWSER_wat" Infos.dat | awk '{ print $2 }'`
   nwat=$(($wat+$wat+$wat))
fi

#
# Overlapping the 100 uncorrelated configurations to generate the ASEC
#
cd MD_ASEC
mv Selected_100.gro Selected_100_full.gro
att=`head -n2 Ordered_100.gro | tail -n1 | awk '{ print $1 }'`
att1=$(($att+1))
att3=$(($att+3))
vol=`head -n $(($att+3)) Ordered_100.gro | tail -n1 | awk '{ print $0 }'`
#   for i in {1..100}
for i in $(eval echo "{1..$confs}")
do
   echo "Generated by MD_ASEC.sh" >> Selected_100.gro
   echo "$shell" >> Selected_100.gro
   head -n $(($att3*$i)) Ordered_100.gro | tail -n $att1 | head -n $(($shell-$nions)) >> Selected_100.gro
   head -n $(($att3*$i)) Ordered_100.gro | tail -n $(($nions+1)) | head -n $nions >> Selected_100.gro
   echo "$vol" >> Selected_100.gro
done

#
# "shell" is the number ot atoms defining the 22 A shell
#
numatoms=$shell
protatoms=`grep "protatoms" ../Infos.dat | awk '{ print $2 }'`

cp $templatedir/ASEC/MD_ASEC.f .
cp ../template_gro2tk .
cp ../template_tk2gro .
echo "Yoe $numatoms $(($protatoms-$nwat))"
sed -i "s|numero|$numatoms|g" MD_ASEC.f
sed -i "s|protatoms|$(($protatoms-$nwat))|g" MD_ASEC.f
sed -i "s|confs|$confs|g" MD_ASEC.f
gfortran MD_ASEC.f -o MD_ASEC.x
./MD_ASEC.x

best=`tail -n1 rms | awk '{ print $1 }'`
#   best=100
head -n $(($att3*$best)) Selected_100_full.gro | tail -n $att3 > Best_Config_full.gro

cd ..

#
# Message to the user
#
cp $templatedir/ASEC/MD_2_QMMM.sh .
./update_infos.sh "Next_script" "MD_2_QMMM" Infos.dat

echo "**************************"
echo ""
echo "Continue with MD_2_QMMM.sh"
echo ""
echo "**************************"

