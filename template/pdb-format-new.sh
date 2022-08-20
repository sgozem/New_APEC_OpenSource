#!/bin/bash
#
# Gromacs (all atoms) --> Tinker PDB conversion
#
# The PDB to use must be the 1st argument of the script
# If it is not found, the script terminates, otherwise some info is given to the user
#
retstereo=`grep "RetStereo" ../Infos.dat | awk '{ print $2 }'`
if [ -z $1 ]; then
   echo ""
   echo " No PDB input file specified."
   echo " You have to type ./pdb-format.sh followed by the PDB file name"
   echo " Terminating..." 
   echo ""
   exit 0
fi
pdbfile=$1
if [ ! -f $pdbfile ]; then
   echo ""
   echo " PDB file not existing. Check if all files are in the same folder."
   echo " Terminating..."
   echo ""
   exit 0
else
   echo ""
   echo " $pdbfile will be modified and saved to final-tk.pdb"
   echo " final-tk.pdb is ready for Tinker pdbxyz conversion"
   echo ""
   echo " Now processing..."
   echo ""
fi
if [[ $retstereo == "nAT" ]]; then
   lyslabel="LYN"
else
   lyslabel="LYS"
fi 
#
# Backing up $pdbfile
#
cp $pdbfile bak_$pdbfile
#
# Generation of replacements vim script
#
echo ':%s/1HH3/ H1 /' > replacements.vim
echo ':%s/2HH3/ H2 /' >> replacements.vim
echo ':%s/3HH3/ H3 /' >> replacements.vim
echo ':%s/OC1/O  /' >> replacements.vim
echo ':%s/OC2/O2 /' >> replacements.vim
echo ':g/SOL/s/ATOM  /HETATM/' >> replacements.vim
echo ':%s/SOL/HOH/' >> replacements.vim
echo ':g/HOH/s/OW/O /' >> replacements.vim
echo ':g/HOH/s/HW1/H1 /' >> replacements.vim
echo ':g/HOH/s/HW2/H2 /' >> replacements.vim
echo ':x' >> replacements.vim
vi -es $pdbfile < replacements.vim
rm replacements.vim
#
# Rearranging lysine and retinal
# Three flags for beginning and ending of retinal, and beginning of water
#
flag=0
flog=0
flug=0
#
# IFS controls the field separator. By putting it like this, echo $line keeps formatting...
#
IFS=""
#
# Reading the input pdbfile and placing retinal in the right position
#
startlys=`awk '{ if ( $3 == "N" && $4 == "RET" ) print $2 }' $pdbfile`
stoplys=`awk '{ if ( $3 == "O" && $4 == "RET" ) print $2 }' $pdbfile`
awk '{ if ( $2 >= '"$startlys"' && $2 <= '"$stoplys"' ) sub($4,"LYS"); print $0 }' $pdbfile > final-tk.pdb
grep RET final-tk.pdb | grep " C." > retinalC
grep RET final-tk.pdb | grep "H" > retinalH
sed -i '/RET/d' final-tk.pdb
 sed -i 's/ CF / C10/;s/ CI / C11/;s/ CJ / C12/;s/ CK / C14/;s/ CM / C15/;s/ CO / C16/;s/ CP / C17/;s/ CQ / C18/;s/ CU / C19/;s/ CX / C20/;s/ATOM  /HETATM/; s/A 296/Z   1/' retinalC
 sed -i 's|HC[0-9][0-9]| HR |;s|[0-9]HC[0-9]| HR |;s|[0-9]HC[A-Z]| HR |;s/ATOM  /HETATM/;s/A 296/Z   1/' retinalH
cat retinalC retinalH > retinalo
echo ":$" > placeret.vim
echo "?ATOM" >> placeret.vim
echo ":r retinalo" >> placeret.vim
echo ":x" >> placeret.vim
vim -es final-tk.pdb < placeret.vim
rm retinalH retinalC retinalo placeret.vim
#
#while read line; do
#     mainline=$line
##
##    $resname is needed for retinal check, along with atom
##    retinal-only part ends when the Lysine backbone N is found
##
#     resname=`echo $mainline | awk '{ print $4 }'`
#     atom=`echo $mainline | awk '{ print $3 }'`
##
##    If flag is not 1, retinal just started. 
##    $lysnum stores the Lysine number for future use
##
#     if [[ $resname == RET && $flag != 1 ]]; then
#        flag=1
#        lysnum=`echo $mainline | awk '{ print $5 }'`
#     fi
##
##    Once N is found, retinal is over as controlled by flog
##
#     if [[ $resname == RET && $atom == N && $flog != 1 ]]; then
#        flog=1
#     fi
##
##   $resname HOH and flug zero mark the beginning of water section
##
#     if [[ $resname == HOH && $flug != 1 ]]; then
#        endprot=`echo $mainline | awk '{ print $2 }'`
#        flug=1
##
##   Now the Lysine part of retinal can be called LYS
##   and retinal C and H stored in file retinal (see below) are added before water
##
#        sed -i "s/RET/$lyslabel/" final-tk.pdb
#        sed -i 's/ CF / C10/;s/ CI / C11/;s/ CJ / C12/;s/ CK / C14/;s/ CM / C15/;s/ CO / C16/;s/ CP / C17/;s/ CQ / C18/;s/ CU / C19/;s/ CX / C20/' retinal
#        grep " C." retinal >> final-tk.pdb
#        grep ".H." retinal | sed 's|HC[0-9][0-9]| HR |;s|[0-9]HC[0-9]| HR |;s|[0-9]HC[A-Z]| HR |' >> final-tk.pdb
#     fi
##
##   If retinal was not met or it is already over, normal PDB creation can proceed
##
#     if [[ $flag == 0 || $flog == 1 ]]; then
#        echo $mainline >> final-tk.pdb
#     else
##
##   When flag is 1 and flog is 0 retinal was found, so it must be stored in the separate file retinal
##
#        echo $mainline | sed 's/ATOM  /HETATM/' >> retinal
#     fi
#done < $pdbfile
#
# Generation of lyso.vim script for lysine modification,
# before applying it to final-tk.pdb
#
#echo 
#echo "/HZ1 LYS   $lysnum" > lyso.vim
#echo "/O   LYS   $lysnum" >> lyso.vim
#echo ":put" >> lyso.vim
#echo ":x" >> lyso.vim
#vi -es final-tk.pdb < lyso.vim
#rm lyso.vim
#
# Changing some hydrogen atoms labels
#
sed -i "s/1HD1/HD11/;s/2HD1/HD12/;s/3HD1/HD13/" final-tk.pdb
sed -i "s/1HD2/HD21/;s/2HD2/HD22/;s/3HD2/HD23/" final-tk.pdb
sed -i "s/1HG1/HG11/;s/2HG1/HG12/;s/3HG1/HG13/" final-tk.pdb
sed -i "s/1HG2/HG21/;s/2HG2/HG22/;s/3HG2/HG23/" final-tk.pdb
#
# Ions need to be HETATM
#
sed -i "/CL/s/ATOM  /HETATM/; /NA/s/ATOM  /HETATM/; /ACI/s/ATOM  /HETATM/" final-tk.pdb
#
# Generating the vim script for chain label replacement if multiple chains are detected
#
multchain=`grep "MultChain" ../Infos.dat | awk '{ print $2 }'`
if [[ $multchain != "YES" && $multchain != "NO" ]]; then
   echo " Cannot read multiple chains information from Infos.dat"
   echo " Are there multiple chains in your protein? (y/n)"
   read answer
   echo ""
   if [[ $answer == "y" ]]; then
      multchain=YES
   else
      multchain=NO
   fi
fi
if [[ $multchain == "YES" ]]; then
   IFS=" "
   lastresa=( $( grep "LastRes" ../Infos.dat ) )
   diffchain=( $( grep "DiffChain" ../Infos.dat ) )
   ngap=${#lastresa[@]}-1
   lastgap=$(($ngap+1))
   lastresa[$lastgap]=1000
   i=1
   for lette in {A..E}; do
       catene[$i]=$lette
       i=$(($i+1))
   done
   IFS=""
   for ((i=1;i<$(($lastgap));i=$(($i+1)))); do
       j=$(($i+1))
       awk '{ if ( $6 > '"${lastresa[$i]}"' ) sub(" '"${catene[$i]}"' "," '"${catene[$j]}"' "); print $0 }' final-tk.pdb > new2.pdb
       mv new2.pdb final-tk.pdb
   done
fi
#
# In any case, retinal must be residue 1 in chain Z, for make Metascript.sh work well
#
retchain=`grep RET final-tk.pdb | head -n 1 | awk '{ print $5 }'`
retnum=`grep RET final-tk.pdb | head -n 1 | awk '{ print $6 }'`
sed -i "s/RET $retchain $retnum/RET Z   1/" final-tk.pdb
#
# Temporary files are removed. Comment these 2 lines for debug
#
cp bak_$pdbfile $pdbfile 
rm bak_$pdbfile
#
# Normal termination
# 
echo " Normal termination - final-tk.pdb was created successfully"
echo " Have a nice day!"
echo ""

