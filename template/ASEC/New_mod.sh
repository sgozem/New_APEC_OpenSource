#!/bin/bash
#
# Copying generic Infos.dat from installation directory
# If Infos.dat already exists, the script terminates to avoid overwriting
#
if [ -f Infos.dat ]; then
   echo ""
   echo " Fatal error - Infos.dat is existing"
   echo " The new project must start in a folder with no other projects"
   echo ""
   echo " Terminating"
   exit 0
fi

templatedir=camino
cp $templatedir/Infos.dat .
Project=Projecto
lysnum=lisina

# Checking if Infos.dat was copied
# If so, it copies the right melacu.prm reading the Tinker version
#
if [ -z Infos.dat ]; then
   echo " Fatal error! - Infos.dat not found"
   echo " Please check the right path of the ARM installation directory in the "fix_pdb.sh" script"
   echo ""
   echo " Terminating..."
   echo ""
   exit 0
else
   echo "Project $Project" >> Infos.dat
fi


while [ -d Step_0 ]; do
      echo " Folder Step_0 found!"
      exit 0
      echo ""
done
mkdir Step_0

# Finding all the PDB files in the current folder
i=1
ext=pdb
for f in *.$ext; do
    filelist[$i]=$f
    i=$(($i+1))
done
i=$(($i-1))

# Variable len counts how many PDB files are in the folder, len+1 is required for ending
# the while loop smoothly. Here the user is asked to select which PDB file has to be used
#
len=${#filelist[*]}
echo ""
echo " Select within $len files $ext:"
echo ""
len=$(($len+1))
i=1
while [ $i -lt $len ]; do
      echo " $i: ${filelist[$i]}"
      let i++
done
echo ""
read choice
if [ -z ${filelist[$choice]} ]; then
   echo " Option unavailable. The program is closing..."
   rm -r $Project
   exit 0
else
   pdbfile=${filelist[$choice]}
   pdb=$(basename $pdbfile .pdb)
   echo " You just selected ${filelist[$choice]}"
   echo ""
fi
cp $pdb.pdb Step_0/${Project}.pdb
if [ -f *.prm ]; then
   cp *.prm Step_0/
fi
mv Infos.dat Step_0/
if [[ -f seqmut ]]; then
   cp Step_0/${Project}.pdb Step_0/wt_${Project}.pdb
   cp seqmut Step_0/
   cp $templatedir/mutate.sh Step_0/
fi

cd Step_0/

# Counts and tells the user the total charge of a system
#
plusca=`awk '{ print $4 " " $5 " " $6 }' ${Project}.pdb | uniq | awk '{ if ( $1 == "LYS" || $1 == "RET" || $1 == "HIP" || $1 == "ARG" ) print $1 }' | wc -l`
pluscb=`awk '{ print $4 " " $5 }' ${Project}.pdb | uniq | awk '{ if ( $1 == "NA" ) print $1 }' | wc -l`
plusc=$(($plusca+$pluscb))
minusca=`awk '{ print $4 " " $5 " " $6 }' ${Project}.pdb | uniq | awk '{ if ( $1 == "ASP" || $1 == "GLU" ) print $1 }' | wc -l`
minuscb=`awk '{ print $4 " " $5 }' ${Project}.pdb | uniq | awk '{ if ( $1 == "CL" || $1 == "ACI" ) print $1 }' | wc -l`
minusc=$(($minusca+$minuscb))
totcarica=$(($plusc-$minusc))

echo "Init_Charge $totcarica" >> Infos.dat

if [[ $totcarica != 0 ]]; then
   echo " The system total charge is $totcarica, which is different than zero"
   echo " Please be sure your PDB file is right"
   echo " Going on..."
   echo ""
else
   echo " The system total charge is $totcarica, everything seems ok"
   echo ""
fi

echo "LysNum $lysnum" >> Infos.dat

# Asking the retinal stereochemistry
#

echo " Please type the retinal type:"
echo " AT  --> all-trans, charge +1"
echo " 11C --> 11-cis, charge +1"
echo " 13C --> 13-cis, charge +1"
echo " nAT --> all-trans, charge 0"
echo " new --> new RESP charges for this model"
read RetStereo
echo ""
case $RetStereo in
     AT) echo "RetStereo AT" >> Infos.dat
     ;;
     11C) echo "RetStereo 11C" >> Infos.dat
     ;;
     13C) echo "RetStereo 13C" >> Infos.dat
     ;;
     nAT) echo "RetStereo nAT" >> Infos.dat
     ;;
     new) echo "RetStereo new" >> Infos.dat
          if [[ -f ../new_charges ]]; then
             res=b
             while [[ $res != "y" && $res != "n" ]]; do
                echo ""
                echo " new_charges file found."
                echo " Will you use it? (y/n)"
                echo ""
                read res
             done
             if [[ $res == "y" ]]; then
                cp ../new_charges .
             else
                echo ""
                echo " Aborting ..."
                echo ""
                exit 0                
             fi
          else
             echo ""
             echo ""
             echo " The RESP charges file (new_charges) is not found. Please provide"
             echo " the full path not including the name and omitting the last /"
             echo ""
             read path
             if [[ -f $path/new_charges ]]; then
                echo ""
                echo " charges found"
                echo ""
                cp $path/new_charges .
             else
                echo ""
                echo " charges not found, Aborting ..."
                echo ""
                exit 0
             fi
          fi
     ;;
     *) echo " Wrong retinal configuration! Aborting..."
        echo ""
        cd ..
        rm -r $Project
        exit 0
esac

# Checking for all the available parameters file and asking for the one to use
#
version=`grep "Tk" Infos.dat | awk '{ print $3 }'`
case $version in
     4.2) cp $templatedir/melacu.prm .
     ;;
     5.1) cp $templatedir/melacu51.prm .
     ;;
     *) echo " Wrong Tinker version detected!"
        echo " New.sh is terminating"
        echo ""
        cd ..
        rm -r $Project
        exit 0
esac
i=1
for prmfile in *.prm; do
    prmlist[$i]=$prmfile
    echo " $i: $prmfile available"
    i=$(($i+1))
done
echo " Please type the number corresponding to the correct prm file"
read choiceprm
if [ -z ${prmlist[$choiceprm]} ]; then
   echo ""
   echo " Wrong prm file selected. Exiting..."
   echo ""
   cd ..
   rm -r $Project
   exit 0
else
   prm=$(basename ${prmlist[$choiceprm]} .prm)
fi
#
# Updating the Infos.dat with prm file name, and CurrCalc field for calculation tracking
#
echo "Parameters $prm" >> Infos.dat
echo "CurrCalc Start" >> Infos.dat

# A folder named as the project is created, and the required files are copied inside
#
#cp $templatedir/NewStep.sh .
#cp $templatedir/smooth_restart.sh .
#cp $templatedir/update_infos.sh .

# Checking for multiple chains in the PDB file. If so, it writes it in Infos.dat 
# for future use
#
grep 'ATOM' ${Project}.pdb > atoms
nchain=`grep 'TER' ${Project}.pdb | wc -l`
chain=A
if [[ $nchain -eq 1 ]]; then
   multchain=no
else
   multchain=yes
   numends=( $( grep -B1 'TER' ${Project}.pdb | grep 'ATOM' | awk '{ print $6 }' ) )
   numstarts=( $( grep -A1 'TER' ${Project}.pdb | grep 'ATOM' | awk '{ print $6 }' ) )
   ngaps=$(($nchain-1))
   for ((i=0;i<$ngaps;i=$(($i+1)))); do
       lastres[$i]=${numends[$i]}
       diffchain[$i]=$((${numstarts[$i]}-${numends[$i]}-1))
   done
#   multchain=yes
#   flag=0
#   j=0
#   while read line && [ $flag -ne $nchain ]; do
#      contr=$chain
#      chain=`echo $line | awk '{ print $5 }'`
#      if [[ $chain != $contr ]]; then
#         flag=$(($flag+1))
#         firstresb=`echo $line | awk '{ print $6 }'`
#         diffchain[$j]=$(($firstresb-$lastresa-1))
#         lastres[$j]=$lastresa
#         j=$(($j+1))
#      fi
#      lastresa=`echo $line | awk '{ print $6 }'`
#   done < atoms
#   rm atoms
fi

# $chain is the chain which retinal belongs to
#
retchain=`grep 'RET' ${Project}.pdb | head -n 1 | awk '{ print $5 }'`
echo "RetChain $retchain" >> Infos.dat
#
# Writing chain info in Infos.dat and running the appropriate pdbxyz command
#
if [[ $multchain == "yes" ]]; then
   echo "MultChain YES" >> Infos.dat
   echo "LastRes ${lastres[@]}" >> Infos.dat
   echo "DiffChain ${diffchain[@]}" >> Infos.dat
else
   echo "MultChain NO" >> Infos.dat
fi
#
# Checking if a list of residues belonging to the retinal cavity is provided
#
#if [[ -f ../cavity ]]; then
#   answer=b
#   while [[ $answer != y && $answer != n ]]; do
#         echo " The cavity file has been found! Do you want to use it? (y/n)"
#         read answer
#         echo ""
#   done
#   if [[ $answer == y ]]; then
#      echo "CavityFile YES" >> Infos.dat
#      cp ../cavity .
#   else
#      echo "CavityFile NO" >> Infos.dat
#   fi
#else
#   echo "CavityFile NO" >> Infos.dat
#fi

#
# Checking the first residue number and writing it in Infos.dat
#
startres=`grep 'ATOM' $Project.pdb | head -n 1 | awk '{ print $6 }'`
echo "StartRes $startres" >> Infos.dat
#
# Retrieving non-standard ionization states for future use
#
for restrano in 'ASH' 'GLH' 'HIE' 'HID' 'LYN'; do
    strange=( $( grep " $restrano " $Project.pdb | awk '{ print $6 }' | uniq ) )
    if [[ ${strange[0]} != '' ]]; then
       echo "$restrano ${strange[@]}" >> Infos.dat
    fi
done
#
# Generate wild type PIR sequence by using Modeller
#
modelchk=`grep "Modeller" Infos.dat | awk '{ print $2 }'`
if [[ $modelchk == YES ]]; then
   cp $templatedir/getpir.py .
   cp $templatedir/mutate_model.py .
   sed -i "s/PROGETTO/$Project/g" getpir.py
   python getpir.py
   echo ""
fi
#
# If seqmut is found, the mutation routine is started, given that either Modeller or SCWRL4 is available
# Such information is retrieved from Infos.dat
#
if [[ -f seqmut ]]; then
   scwrlchk=`grep "SCWRL4" Infos.dat | awk '{ print $2 }'`
   modelchk=`grep "Modeller" Infos.dat | awk '{ print $2 }'`
   if [[ $scwrlchk == YES || $modelchk == YES ]]; then 
      ./mutate.sh $Project $scwrlchk $modelchk
   else
      echo " You don't have any software to perform mutations!"
      echo " Please delete seqmut and try again"
      echo " Aborting..."
      echo ""
      exit 0
   fi
   errid=`grep 'mutate_errid' arm.err | awk '{ print $2 }'`
   if [[ $errid -eq 0 ]]; then
      echo " The following mutations has been inserted:"
      cat seqmut
      echo ""
   else
      echo " Problem in mutate.sh, no mutations will be performed."
      echo " Please type y to go on with the wild type, or n to abort"
      read scelta
      echo ""
      case $scelta in
           y) echo " Going on with the wild type..."
              mv wt_${Project}.pdb $Project.pdb
           ;;
           n) echo " Aborting..."
              echo ""
              exit 0
           ;;
      esac
   fi
else
   echo " No mutants requested, going on with the wild type.."
   echo ""
fi
#
# Messages to the user
#

cp $templatedir/ASEC/Init_ASEC.sh .

echo " Folder Step_0 created! cd Step_0/, then ./Init_ASEC.sh"
echo ""


