#!/bin/bash
#
# Reading information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
solvent=`grep "SolventBox" Infos.dat | awk '{ print $2 }'`
tmpdir=`grep "tempdir" Infos.dat | awk '{ print $2 }'`

#
# Collecting data to run the QM/MM calculations
#
cp MD_ASEC/list_tk.dat calculations/${Project}_OptSCF
cp MD_ASEC/ASEC_tk.xyz calculations/${Project}_OptSCF
cp $templatedir/ASEC/ASEC.f calculations/${Project}_OptSCF
cp $prm.prm calculations/${Project}_OptSCF
cp $templatedir/ASEC/New_parameters.f calculations/${Project}_OptSCF

cd calculations/${Project}_OptSCF
cp ${Project}_OptSCF.xyz ${Project}_OptSCF_old.xyz
mv ${Project}_OptSCF.xyz coordinates_tk.xyz

#
# Generating the final ASEC superconfigurations including the
# atom type of the ASEC pseudo-atoms
#
if [[ -f ../../chargefxx ]]; then
   cp ../../chargefxx .
   fxx=`grep -c "CHARGE" chargefxx | awk '{ print $1 }'`
   sed -i "s|tailall|$fxx|g" ASEC.f
   sed -i "s|tailall|$fxx|g" New_parameters.f
else
   sed -i "s|tailall|0|g" ASEC.f
   sed -i "s|tailall|0|g" New_parameters.f
fi

shell=`grep "Shell" ../../Infos.dat | awk '{ print $2 }'`

#
# $shell+1 for the link atom
#
sed -i "s|numero|$(($shell+1))|g" ASEC.f
gfortran ASEC.f -o ASEC.x
./ASEC.x
mv new_coordinates_tk.xyz ${Project}_OptSCF.xyz

#
# This section is for obtaining the new force fiel parameters
# of the ASEC pseudo-atoms
#

grep "atom     " $prm.prm > atom.dat
numindx2=`grep -c "atom     " atom.dat`
sed -i "s|atomos|$numindx2|g" New_parameters.f

grep "charge   " $prm.prm > charges.dat
numcharges=`grep -c "charge   " charges.dat`
sed -i "s|cargas|$numcharges|g" New_parameters.f

#
# $shell+1 for the link atom
#
sed -i "s|numero|$(($shell+1))|g" New_parameters.f

grep "vdw      " $prm.prm > vdw.dat
numvdw=`grep -c "vdw      " vdw.dat`
sed -i "s|vander|$numvdw|g" New_parameters.f
sed -i "s|chnorm=40|chnorm=100|g" New_parameters.f

gfortran New_parameters.f -o New_parameters.x
./New_parameters.x

atom=`grep -n "atom     " $prm.prm | awk '{print $1}' FS=":" | tail -n1`
head -n$atom $prm.prm > temp1
cat temp1 new_atom.dat > temp2

vdw=`grep -n "vdw      " $prm.prm | awk '{print $1}' FS=":" | tail -n1`
head -n$vdw $prm.prm | tail -n$(($vdw-$atom)) >> temp2
cat temp2 new_vdw.dat > temp3

charge=`grep -n "charge   " $prm.prm | awk '{print $1}' FS=":" | tail -n1`
head -n$charge $prm.prm | tail -n$(($charge-$vdw)) >> temp3
cat temp3 new_charges.dat > temp4

last=`wc -l $prm.prm | awk '{print $1}'`
tail -n$(($last-$charge)) $prm.prm >> temp4

mv $prm.prm old_$prm.prm 
mv temp4 $prm.prm

#rm ASEC_tk.xyz list_tk.dat coordinates_tk.xyz ASEC.x ASEC.f
rm temp1 temp2 temp3 new_atom.dat new_charges.dat new_vdw.dat atom.dat charges.dat vdw.dat New_parameters.x New_parameters.f

#
# Submiting the job
#
#TMPFILE=`mktemp -d /scratch/photon_XXXXXX`
#../../update_infos.sh "tempdir" $TMPFILE ../../Infos.dat
#sed -i "s|TEMPFOLDER|$TMPFILE|g" molcas-job.sh
#cp -r * $TMPFILE
#current=$PWD
#cd $TMPFILE
sbatch molcas-job.sh
#cd $current

#
# Message to the user
#
cd ..
cp $templatedir/ASEC/Molcami2_mod.sh .
../update_infos.sh "Next_script" "Molcami2_mod.sh" ../Infos.dat
echo ""
echo "******************************************************"
echo ""
echo " After complete the OptSCF/ANO-L-VDZ, run Molcami2_mod.sh"
echo ""
echo "******************************************************"


