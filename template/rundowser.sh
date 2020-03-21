#!/bin/bash
#
#
# Script for running Dowser - called by NewStep.sh
# Usage: ./rundowser.sh [ project name ] [ Tinker folder ] [ prm filename ]
#
Project=$1
tinkerdir=$2
prm=$3
pathdow=`grep "DowPath" ../../Infos.dat | awk '{ print $2 }'`
multi=`grep "MultChain" ../../Infos.dat | awk '{ print $2 }'`
startres=`grep "StartRes" ../../Infos.dat | awk '{ print $2 }'`
#
# Grepping water from PDB. Alanines are required to make Tinker recognize water
#
grep 'ALA' $Project.pdb > H2O.pdb
grep 'HOH' $Project.pdb >> H2O.pdb
#
# Adding hydrogen atoms with Tinker by converting the PDB into XYZ and vice versa
#

echo "Yoe Debug 1"

if [[ $multi == "YES" ]]; then
$tinkerdir/pdbxyz H2O.pdb <<EOF
ALL
${prm}.prm
EOF
else
$tinkerdir/pdbxyz H2O.pdb <<EOF
${prm}.prm
EOF
fi
$tinkerdir/xyzpdb H2O.xyz <<EOF
${prm}.prm
EOF

echo "Yoe Debug 2"

#
# Removing the CONECT tags and renaming WAT as HOH
#
echo -e 'g/CONECT/d' >> aggiusta
echo -e ':%s/WAT/HOH/g' >> aggiusta
echo -e ':x' >> aggiusta
vi -es  H2O.pdb_2 < aggiusta
#
# PdbFormatter renames the water hydrogen atoms, from HT to H1 and H2
# output file is file.out
#
./PdbFormatter.py
rm PdbFormatter.py
#
# Removing water from the protein PDB
#
echo -e 'g/HOH/d' >> toglih2o
echo -e 'g/END/d' >> toglih2o
echo -e ':x' >> toglih2o
vi -es  $Project.pdb < toglih2o
#
# Adding file.out and END to the protein PDB
# file.out is the PdbFormattery.py output
#
cat file.out >> $Project.pdb
echo 'END' >> $Project.pdb
echo -e ':%s/ OT / OW /g' >> aggiusta2
echo -e ':x' >> aggiusta2
vi -es  $Project.pdb < aggiusta2
#
# Cleaning up scratch files from previous steps
#
#  Yoe rm toglih2o aggiusta* H2O.* file.out
#
# Executing pdb-to-dow to reformat the PDB in a suitable way for Dowser 
#
./pdb-to-dow.sh $Project.pdb
if [[ -f ../../arm.err ]]; then
   checkpdbdow=`grep 'pdb-to-dow' ../../arm.err | awk '{ print $2 }'`
fi
if [[ $checkpdbdow -ne 0 ]]; then
   echo " An error occurred in pdb-to-dow.sh. I cannot go on..."
   echo ""
   echo "rundowser.sh 1 PDBtoDowProblem" >> ../../arm.err
   exit 0
fi
#
# Asking the user about the use of crystal waters in Dowser guess
#

echo ""
echo ""
echo ""
echo "####################################################################################"
echo ""
echo " Dowser is a computational code used in this protocole to deal with the water  "
echo " molecules and cavities fund in the PDB crystallographic structure." 
echo " There are two available options: "
echo ""
echo " 1 Consider just the water molecules found in the PDB crystallographic structure"
echo "   removing the non-bounded one."
echo " 2 Desconsider all water molecules fund in the PDB crystallographic structure and"
echo "   add water molecules into the empty cavities based on an energetic and steric "
echo "   criterium"
echo ""
echo "  In this case option 1 will be used"
echo "####################################################################################"
echo ""
echo ""
sleep 10

answer="y"
#while [[ $answer != y && $answer != n ]]; do
#      echo " Type y to use only crystal waters or type n to add waters by Dowser?"
#      read answer
#done
echo ""
if [[ $answer == y ]]; then
   argo="-onlyxtalwater"
else
   argo=""
fi
#
# Dowser execution, $pathdow includes the full dowser path from Infos.dat
#
$pathdow dowser.pdb $argo
#
# reform.pdb is the Dowser output, which must be adapted for the subsequent conversion to Gromacs
#
cp reform.pdb reform.new.pdb
#
# Renaming the retinal labels according to Gromacs format
#
sed -i 's/C10/CF /' reform.new.pdb
sed -i 's/C11/CI /' reform.new.pdb
sed -i 's/C12/CJ /' reform.new.pdb
sed -i 's/C14/CK /' reform.new.pdb
sed -i 's/C15/CM /' reform.new.pdb
sed -i 's/C16/CO /' reform.new.pdb
sed -i 's/C17/CP /' reform.new.pdb
sed -i 's/C18/CQ /' reform.new.pdb
sed -i 's/C19/CU /' reform.new.pdb
sed -i 's/C20/CX /' reform.new.pdb
#
# At the end of the PDB, all the water from Dowser and chloride ions are added
#
echo TER >> reform.new.pdb
cat dowserwat_all.pdb >> reform.new.pdb
cat CL >> reform.new.pdb
cat NA >> reform.new.pdb
cat ACI >> reform.new.pdb
#
# Adding the correct chain labels to all residues
#
echo -e ':%s/ATOM................./&A/g' >> mer
echo -e ':%s/HOH W/HOH A/g' >> mer
echo -e '%s/ A  / A /g' >> mer
echo -e ':x' >> mer
vi -es reform.new.pdb < mer
rm mer
#
# Backing up reform.new.pdb and execution of H renaming script
#
cp reform.new.pdb reform.new.pdb.Old
./yesH-tk-to-gro.sh reform.new
#
# Reshifting the numbering in the PDB when multiple chains are present
#
sed -i "/REMARK/d" reform.new.pdb
if [[ $multi == "YES" ]]; then
   lastresa=( $( grep "LastRes" ../../Infos.dat ) )
   diffchain=( $( grep "DiffChain" ../../Infos.dat ) )
   ngap=${#lastresa[@]}-1
   for ((j=1;j<=$(($ngap));j=$(($j+1)))); do
       shifter=$((${diffchain[$j]}))
       awk '{ if ( $6 > '"${lastresa[$j]}"' ) { newnum = $6 + '"$shifter"'}  else { newnum = $6}; printf("%3s\n", newnum) }' reform.new.pdb > resnum
       cut -c1-22 reform.new.pdb > testa
       cut -c28-106 reform.new.pdb > coda
       paste -d" " testa resnum coda > numok_reform.new.pdb
       mv numok_reform.new.pdb reform.new.pdb
       rm testa resnum coda
   done
   mv reform.new.pdb nuovissimo.pdb
#   lastresa=( $( grep "LastRes" ../../Infos.dat ) )
#   diffchain=( $( grep "DiffChain" ../../Infos.dat ) )
#   ngap=${#lastresa[@]}-1
#   j=1
#   diffs=$(($startres-1))
#   lastres=${lastresa[$j]}
#   k=2
#   IFS=""
#   while read line; do
#         resnum=`echo $line | awk '{ print $6 }'`
#         restype=`echo $line | awk '{ print $4 }'`
#         vero=$(($lastres-$diffs))
#               echo $resnum $diffs $vero >> cazzone
#         if [[ $resnum -gt $vero ]]; then
#            diffs=$(($diffs+${diffchain[$j]}))
#            if [[ $k -le $ngap ]]; then
#               lastres=${lastresa[$k]}
#               j=$k
#               k=$(($k+1))
#            else
#               lastres=1000
#            fi
#         fi
#         resnum1=$resnum
#         resnum=$(($resnum1+$diffs))
#         if [[ $restype != "HOH" ]]; then
#            testo=`echo $line | awk '{ if ( $2 == $6 ) print "yes" }'`
#            if [[ $testo != "yes" ]]; then
#               if [[ $resnum1 -lt 100 && $resnum -ge 100 ]] || [[ $resnum1 -lt 10 && $resnum -ge 10 ]]; then
#                  echo $line | sed "s/ $resnum1 /$resnum /" >> nuovissimo.pdb
#               else
#                  echo $line | sed "s/ $resnum1 / $resnum /" >> nuovissimo.pdb
#               fi
#            else
#                  echo $line | sed "s/ $resnum1 /$resnum /2" >> nuovissimo.pdb
#            fi
#         else
#            echo $line >> nuovissimo.pdb
#         fi
#   done < reform.new.pdb
else
  mv reform.new.pdb nuovissimo.pdb
fi
#
# The PDB is ready to undergo Gromacs conversion
#
cp nuovissimo.pdb $Project.pdb
#
# If the initial residue number is different from 1, Dowser put it as 1
# The following section fixes residue numbering consistently with the initial PDB
#
firstnum=`grep '^ATOM' ../$Project.pdb | head -n 1 | awk '{ print $6 }'`
if [[ $firstnum -gt 1 ]]; then
   shifter=$(($firstnum-1))
   awk '{ if ( $1 != HETATM ) { newnum = $6 + '"$shifter"'}  else { newnum = $6}; printf("%3s\n", newnum) }' $Project.pdb > resnum
   cut -c1-22 $Project.pdb > testa
   cut -c28-106 $Project.pdb > coda
   paste -d" " testa resnum coda > numok_$Project.pdb
   mv $Project.pdb dowserfinal_$Project.pdb
   mv numok_$Project.pdb $Project.pdb
   rm testa resnum coda
fi
echo "rundowser.sh 0 OK" >> ../../arm.err
