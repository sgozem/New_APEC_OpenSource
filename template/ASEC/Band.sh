#!/bin/bash
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`

#cd calculations
if [ -d Band_width ]; then
   echo ""
   echo " Folder calculations already exists. Please check it and remove if necessary"
   echo " Terminating..."
   echo ""
   exit 0
fi

confs=100

mkdir Band_width
cd Band_width

cp ../MD_ASEC/Selected_100.gro .
cp ../$prm.prm .

numat=`head -n2 Selected_100.gro | tail -n1 | awk '{ print $1 }'`

cp ../calculations/${Project}_6-31G_Opt/${Project}_6-31G_Opt_old.xyz .
numatt=`head -n1 ${Project}_6-31G_Opt_old.xyz | awk '{ print $1 }'`
head -n $(($numatt+1)) ${Project}_6-31G_Opt_old.xyz | tail -n1 > la

cav=b
while [[ $cav != "y" && $cav != "n" ]]; do
   echo ""
   echo " In addition to the structures for computing the Absorption band,"
   echo " you can define a cavity for each of the generated structures"
   echo " for a subsequent Molcas/Tinker Dynamics. Do you want to do it? (y/n)"
   echo ""
   read cav
done

right=b
if [[ $cav == "y" ]]; then
   while [[ $right != "y" && $right != "n" ]]; do
      echo ""
      echo " What is the size of the cavity in Angstroms"
      echo ""
      read size
      echo ""
      echo " Relax also the backbone? (y/n)"
      echo ""
      read backbone
      
      echo ""
      echo " So, $size Angstrom cavity with relax backbone \"$backbone\". Right? (y/n)"
      echo ""
      read right
   done
   mkdir ../Tinker_${size}A_cavity

   echo ""
   echo ""
   echo " The folder \"Tinker_${size}A_cavity\" will be created with the $confs xyz configurations"
   echo " and the corresponding key files"
   echo ""
   sleep 3

   cp $prm.prm ../Tinker_${size}_A_cavity
   if [[ $backbone == "y" ]]; then
       selection="same residue as (all within $size of (resname RET and not name N H CA HA C O)) and not (resname RET)"
   else
      selection="(same residue as (all within $size of (resname RET and not name N H CA HA C O))) and (sidechain or water or (resname ACI ACH)) and not resname RET"
   fi
fi

for i in $(eval echo "{1..$confs}"); do

   # to convert gro into pdb into Tinker xyz
   #   
   head -n $(($i*($numat+3))) Selected_100.gro | tail -n $(($numat+3)) > final-$Project.gro

 if [[ $cav == "y" ]]; then
   cp ../template_gro2tk .
   # TCL script for VMD: open file, apply selection, save the serial numbers into a file
   #
   echo -e "mol new final-$Project.gro type gro" > ndxsel.tcl
   echo -e "mol delrep 0 top" >> ndxsel.tcl
   riga1="set babbeo [ atomselect top \"$selection\" ]"
   echo -e "$riga1" >> ndxsel.tcl
   echo -e 'set noah [$babbeo get serial]' >> ndxsel.tcl
   riga3="set filename grodyn.dat"
   echo -e "$riga3" >> ndxsel.tcl
   echo -e 'set fileId [open $filename "w"]' >> ndxsel.tcl
   echo -e 'puts -nonewline $fileId $noah' >> ndxsel.tcl
   echo -e 'close $fileId' >> ndxsel.tcl
   echo -e "exit" >> ndxsel.tcl
   vmd -e ndxsel.tcl -dispdev text

   cp grodyn.dat list

   if [[ -f list_tk ]]; then
      rm list_tk
   fi
   num=`awk '{print NF}' list`
   for k in $(eval echo "{1..$num}")
   do
     atomgro=`awk -v j=$k '{ print $j }' list`
     atomtk=`head -n $(($atomgro+1)) template_gro2tk | tail -n1 | awk '{ print $2 }'`
     echo "ACTIVE $atomtk" >> list_tk
   done
#   rm list
 fi  

   cp $templatedir/pdb-format-new.sh .
   $gropath/gmx editconf -f final-${Project}.gro -o final-$Project.pdb -label A
   ./pdb-format-new.sh final-$Project.pdb
   # pdbxyz conversion
   #
   mv final-tk.pdb $Project-tk.pdb
   $tinkerdir/pdbxyz $Project-tk.pdb << EOF
ALL
../$prm
EOF
   cat $Project-tk.xyz la > a
   mv a ${Project}_conf_${i}.xyz
   line=`head -n1 ${Project}_conf_${i}.xyz | awk '{ print $0 }'`

   sed -i "s/$line/$numatt/" ${Project}_conf_${i}.xyz

   cp ../calculations/CASPT2_ipea_0/${Project}_CASPT2_0.key ${Project}_conf_${i}.key
   if [[ $cav == "y" ]]; then
      cp ${Project}_conf_${i}.xyz ../Tinker_${size}A_cavity
      cat ${Project}_conf_${i}.key list_tk > a
      mv a ../Tinker_${size}A_cavity/${Project}_conf_${i}.key
      cp list_tk ../Tinker_${size}A_cavity
      cd ../Tinker_${size}A_cavity
         sed -i "s/QMMM 63/QMMM $((63+$num))/g" ${Project}_conf_${i}.key
         line1=`grep -n "QMMM-ELECTROSTATICS" ${Project}_conf_${i}.key | cut -d : -f 1`
         head -n $(($line1-1)) ${Project}_conf_${i}.key > headd
         sed -i "s/ACTIVE/MM/g" list_tk
         cat headd list_tk > a
         mv a headd
         todo=`wc -l ${Project}_conf_${i}.key | awk '{ print $1 }'`
         tail -n $(($todo-$line1-1)) ${Project}_conf_${i}.key > bottom
         cat headd bottom > a
         mv a ${Project}_conf_${i}.key
      cd ../Band_width
   fi
   cp ../calculations/CASPT2_ipea_0/${Project}_CASPT2_0.input ${Project}_conf_${i}.input
   rm $Project-tk.pdb $Project-tk.seq final-$Project.pdb $Project-tk.xyz final-$Project.gro

   mod=$(($i % 5))
   if [[ $mod -eq 0 ]]; then
      for j in {4..0}; do   
         cp $templatedir/molcas-job.sh molcas-job_$(($i-$j)).sh

         sed -i "s|NOMEPROGETTO|${Project}_conf_$(($i-$j))|" molcas-job_$(($i-$j)).sh
         sed -i "s|MEMTOT|23000|" molcas-job_$(($i-$j)).sh
         sed -i "s|MEMORIA|20000|" molcas-job_$(($i-$j)).sh
         sed -i "s|hh:00:00|20:00:00|" molcas-job_$(($i-$j)).sh
         sed -i "s|ppn=16|ppn=4|" molcas-job_$(($i-$j)).sh

         sed -i "/#PBS -m ae/d" molcas-job_$(($i-$j)).sh
    
#         qsub molcas-job_$(($i-$j)).sh

         echo ""
         echo " Sumbitted configuration $(($i-$j)) ..."
         echo ""
         echo ""
      done
   fi
done

