#!/bin/bash
#
# 
# Allow the user to select between Modeller and Scwrl4 to perform mutations
#
Project=$1
scwrlchk=$2
modelchk=$3
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
startres=`grep "StartRes" Infos.dat | awk '{ print $2 }'`
multi=`grep "MultChain" Infos.dat | awk '{ print $2 }'`

if [[ $scwrlchk == YES && $modelchk == YES ]]; then
   choice=3
   while [[ $choice -ne 1 ]]  && [[ $choice -ne 2 ]]; do
         echo " Please select the program you want to use to performe the mutations:"
         echo " 1. - Modeller"
         echo " 2. - Scwrl4"
         read choice
         echo ""
   done
else
   if [[ $scwrlchk == YES ]]; then
       echo " SCWRL4 is installed, while Modeller is not. Using SCWRL4..."
       echo ""
       choice=2
   else
       echo " Modeller is installed, while SCWRL4 is not. Using Modeller..."
       echo ""
       choice=1
   fi
fi
#
# According to user's choice, the Modeller python routine or the runscwrl.sh
# launch script are copied into the current folder
#
case $choice in
     1) cp $templatedir/mutate_model.py .
        echo " Using Modeller..."
        echo ""
     ;;
     2) echo " Using SCWRL4..."
        echo ""
     ;;
     *) echo " A problem occurred! Aborting..."
        echo "" 
        exit 0
     ;;
esac
#
# seqmut is the file including the point mutations. Multiple mutations
# must be on different lines. Replacements must be written as:
# [oldres][number][newres]
# where oldres and newres can be 1-letter or 3-letter aminoacid standard codes.
# The names are case insensitive
#
# Examples:
# SER33thr
# K89L
# k143W
# ser99Pro
#
# Mixed 1-letter and 3-letter codes as in ASP56K will not work!!
#
# mutation is a vector of strings with all mutants, converted to the uppercase
# num is the number of its elements, i.e. the number of mutations
#
mutation=( $( cat seqmut | tr "[:lower:]" "[:upper:]" ) )
num=${#mutation[@]}
#
# List the detected substitutions
#
echo " seqmut includes the following $num single point mutations:"
echo " ${mutation[@]}"
echo ""
#
# Loop over each replacement
# Every string is checked for 1- or 3-letters codes and the procedure is 
# modified accordingly. In the case of 1-letter codes, the aminoconv script is copied and called
#
restutti="GLY ALA VAL LEU ILE SER THR CYS MET PRO HIS ARG ASN GLN GLU ASP PHE TRP TYR LYS ASH GLH HID HIE LYN"
for ((i=0;$i<$num;i=$(($i+1)))); do
    length=${#mutation[$i]}
    curmut=${mutation[$i]}
    controllo=`expr match "$curmut" '[A-Z][0-9]'`
    if [[ $controllo -eq 0 ]]; then
       last3=$(($length-3))
       wtres=${curmut:0:3}
       sidemut=${curmut:$last3:3}
       if [[ $sidemut == ASH || $sidemut == GLH || $sidemut == HID || $sidemut == HIE || $sidemut == LYN  ]]; then
          echo $sidemut > NeutralMut
       fi
       reschk=`echo $restutti | grep "$sidemut" | wc -l`
       if [[ $reschk -eq 0 ]]; then
          echo " The residue you want to introduce is not recognized! Please check your mutation!"
          echo " Aborting..."
          echo ""
          echo "mutate.sh 1 NewResidueNotRecognized" >> arm.err
          exit 0
       fi
    else
       last1=$(($length-1))
       wtres1=${curmut:0:1}
       sidemut1=${curmut:$last1:1}
       cp $templatedir/aminoconv.sh .
       ./aminoconv.sh $wtres1 $sidemut1 3
       wtres=`cat wtres`
       sidemut=`cat sidemut`
       rm wtres sidemut
    fi
#
# After the if selection, curmut include the full string, wtres and
# sidemut the 3-letter codes for the native and replacing residues respectively
# the digits in curmut are the residue number, numwt
# if the label of the native residue matches the replacement, the mutation is
# wrong and the script exits
#
    numwt=`echo $curmut | sed 's/[^0-9]//g'`
    checkres=`awk '{ if ( $6 == "'"$numwt"'" ) print $4 }' $Project.pdb | head -n 1`
    if [[ $checkres == $sidemut ]]; then
       echo " The new residue you want to insert is already present in the structure!"
       echo " Do you want to rearrange it? (y/n)"
       read answer
       echo ""
       if [[ $answer == n ]]; then
          echo " No further operation! Aborting..."
          echo "mutate.sh 2 NewResidueIdenticalOld" >> arm.err 
          exit 0
       fi
    else
       if [[ $checkres != $wtres ]]; then
          echo " Wild type residue not found! Please check your mutation!"
          echo " Aborting..."
          echo ""
          echo "mutate.sh 3 WrongWildTypeResidue" >> arm.err
          exit 0
       fi
    fi
#
# After reading the chain label, required by Modeller, it runs the Modeller Python
# script if the user selected Modeller.
# Since Modeller changes the order of backbone atoms, the newly inserted residue
# is taken from its output file, its atom are shifted by shiftatm.vim and it is
# placed back in the PDB input file in place of the native one (replacemut.vim)
#
    chainwt=`awk '{ if ( $6 == "'"$numwt"'" ) print $5 }' $Project.pdb | head -n 1`
    if [[ $choice -eq 1 ]]; then
       python mutate_model.py $Project $numwt $sidemut $chainwt > mut_$Project.log
       awk '{ if ( $6 == '"$numwt"' ) print $0 }' "$Project""_""$sidemut""$numwt"".pdb" > nuovomut.pdb
       echo ":/ C   /" > shiftatm.vim
       echo ':.,$d' >> shiftatm.vim
       echo ":/ CA  /" >> shiftatm.vim
       echo ":put" >> shiftatm.vim
       echo ":x" >> shiftatm.vim
       vim -es nuovomut.pdb < shiftatm.vim
       awk '{ if ( $6 != '"$numwt"' ) print $0 }' $Project.pdb > fetecchia.pdb
       mv fetecchia.pdb $Project.pdb
       numeno=$(($numwt-1))
       echo ":$" > replacemut.vim
       echo "? $numeno " >> replacemut.vim
       echo ":r nuovomut.pdb" >> replacemut.vim
       echo ":x" >> replacemut.vim
       vim -es $Project.pdb < replacemut.vim
       rm replacemut.vim shiftatm.vim nuovomut.pdb "$Project""_""$sidemut""$numwt"".pdb"
    else
#
# Scwrl4 branch - after copying $Project.pdb into scwrl.pdb, all operations are
# performed on such file. Re-labelling of non standard ionization states was already done before
#
       cp $Project.pdb scwrl.pdb
       sed -i 's/ HIE / HIS /g' scwrl.pdb
       sed -i 's/ HID / HIS /g' scwrl.pdb
       sed -i 's/ ASH / ASP /g' scwrl.pdb
       sed -i 's/ GLH / GLU /g' scwrl.pdb
       sed -i 's/ LYN / LYS /g' scwrl.pdb
#
# getpir.py is the Python Modeller script to write residue sequence in PIR format from any PDB
# The obtained wild type sequence will be modified with the mutation and given to SCWRL4 as an input
#
       if [[ $modelchk == YES && $choice -eq 1 ]]; then
          cp $templatedir/getpir.py .
          sed -i "s/PROGETTO/scwrl/g" getpir.py
          python getpir.py
#
# Some manipulation on PIR files: first 3 rows deleted, all residues put on 1 line only, removed the
# '*' from the end...basically it is not a PIR file anymore, except for the '/' which indicates gaps
#
          echo ":1,3d" > cleanpir.vim
          echo ":%s/\n//" >> cleanpir.vim
          echo "s/\*//" >> cleanpir.vim
          echo ":x" >> cleanpir.vim
          vim -es scwrl.pir < cleanpir.vim
          rm cleanpir.vim
       else   # this is needed for scwrl4
          cp $templatedir/aminoconv.sh .
          cp $templatedir/getpir.sh .
          sed -i "s/PROGETTO/scwrl/g" getpir.sh     
          ./getpir.sh 
       fi
       grep 'HETATM' scwrl.pdb > retwat.pdb 
#
# SCWRL4 treats all lowercase residues as frozen, so the original uppercase sequence must be made lowercase
#
       cat scwrl.pir | tr "[:upper:]" "[:lower:]" > okscwrl.pir
       sequel=`cat okscwrl.pir`
#
# From Infos.dat the number of gaps is counted, so that the array of residue symbols has the correct sequence
# numbering, for example if there is a gap of 11 residues after ASP 54, the following residue will be 66,
# so the array sequenza must have no symbols in positions 55, 56 and so on
#
# Bug fixed by Yoe
#
       if [[ $multi == "YES" ]]; then
          gap=( $( grep DiffChain Infos.dat ) )
       fi
       if [[ $multi == "NO" ]] && [[ $startres -ne 1 ]]; then
          gap[1]=$(($startres-1))
       fi
       quanti=${#sequel}
       k=1
       p=1
       for ((j=0;$j<$quanti;j=$(($j+1)))); do
           if [[ ${sequel:$j:1} != "/" ]]; then
              sequenza[$k]=`echo ${sequel:$j:1}`
           else
              k=$(($j+${gap[$p]}))
              quanti=$(($quanti+${gap[$p]}))
              p=$(($p+1))
           fi
           k=$(($k+1))
       done
#
# aminoconv.sh converts 3-letter codes to 1-letter and viceversa 
# Regardless the notation used in seqmut, the uppercase 3-letter codes are used,
# therefore we need to convert them for SCWRL4
#
       cp $templatedir/aminoconv.sh .
       ./aminoconv.sh $wtres $sidemut 1 
       wtres1=`cat wtres | tr "[:upper:]" "[:lower:]"`
       sidemut1=`cat sidemut`
       rm wtres sidemut
#
# Mutation is performed on the residue array sequenza, by changing the 1-letter code
# corresponding to the residue number from seqmut (which is stored in $nuwt)
#
       if [[ ${sequenza[$numwt]} == $wtres1 ]]; then
          sequenza[$numwt]=$sidemut1
       else
          echo " Wild type residue not found! Please check your mutation!"
          echo " Aborting..."
          echo ""
          echo "mutate.sh 3 WrongWildTypeResidue" >> arm.err
          exit 0
       fi
#
# This loop writes the mutated PIR sequence, simply by writing the array content on the same line in mutated.pir
#
       if [[ -f mutated.pir ]]; then
          rm mutated.pir
       fi
       for j in $(seq 1 $quanti); do
          echo -n ${sequenza[$j]} >> mutated.pir
       done 
#
# Call to SCWRL4: pdb, mutated PIR sequence and retwat.pdb (the HETATM file, needed to put vdW constraints)
# the -h flag prevents the hydrogen atoms addition. putret.vim put retinal and water in the new pdb file
# then some cleanup is performed
#
       Scwrl4 -i scwrl.pdb -o "$Project""_""$sidemut""$numwt"".pdb" -s mutated.pir -h -f retwat.pdb > scwrl4_mut.log
       awk '{ if ( $6 == '"$numwt"' ) print $0 }' "$Project""_""$sidemut""$numwt"".pdb" > nuovomut.pdb
       if [[ -f NeutralMut ]]; then
          curside=`awk '{ print $4 }' nuovomut.pdb | uniq`
          truemut=`cat NeutralMut`
          sed -i "/ $numwt /s/$curside/$truemut/" nuovomut.pdb
          rm NeutralMut
       fi
#       numeno=$(($numwt-1))

       awk '{ if ( $6 != '"$numwt"' ) print $0 }' $Project.pdb > fetecchia.pdb
       mv fetecchia.pdb $Project.pdb

#
# Bug fixed by Yoe
#

       awk '{ if ( $6 == '"$numwt"' ) print $0 }' "scwrl.pdb" > temp
       res=`head -n1 temp | awk '{ print $4 }'`
       chain=`head -n1 temp | awk '{ print $5 }'`
       line=`grep -n -m1 "$res $chain $numwt" scwrl.pdb | awk '{ print $1 }' FS=":"`
       prev_res=`head -n $(($line-1)) scwrl.pdb | tail -n1 | awk '{ print $4 }'`
       prev_chain=`head -n $(($line-1)) scwrl.pdb | tail -n1 | awk '{ print $5 }'`

       echo ":$" > replacemut.vim
       echo "? $prev_res $prev_chain $(($numwt-1))" >> replacemut.vim
       echo ":r nuovomut.pdb" >> replacemut.vim
       echo ":x" >> replacemut.vim
       vim -es $Project.pdb < replacemut.vim
       rm replacemut.vim nuovomut.pdb "$Project""_""$sidemut""$numwt"".pdb"
#
# Putting back the non-standard ionization states in the mutated pdb file
#
#       for restrano in 'ASH' 'GLH' 'HIE' 'HID' 'LYN'; do
#           strange=( $( grep "$restrano" Infos.dat ) )
#           numstr=${#strange[@]}
#           j=0
#           while [[ $j -lt $numstr ]]; do
#                 j=$(($j+1))
#                 if [[ ${strange[$j]} -ne $numwt ]]; then
#                    awk '{ if ( $6 == "'"${strange[$j]}"'"  &&  match($0,"TER") == 0 ) sub($4, "'"${strange[0]}"'"); print $0 }' "$Project""_""$sidemut""$numwt"".pdb" > scratch.pdb
#                    mv scratch.pdb "$Project""_""$sidemut""$numwt"".pdb"
#                 fi
#           done
#       done
#       mv "$Project""_""$sidemut""$numwt"".pdb" $Project.pdb
    fi
done
echo "mutate.sh 0 OK" >> arm.err

