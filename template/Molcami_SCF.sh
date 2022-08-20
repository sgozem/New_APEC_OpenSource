#!/bin/bash
#
# Retrieving the needed information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
retstereo=`grep "RetStereo" Infos.dat | awk '{ print $2 }'`
version=`grep "Tk" Infos.dat | awk '{ print $3 }'`
#
# Creating the directory calculations
# If calculations already exists, the script asks if this is a restart
#
echo ""
echo " Using $Project.xyz and $prm.prm"
echo ""
if [ -d calculations ]; then
   ./smooth_restart.sh calculations "Do you want to re-run all the QM/MM calculations? (y/n)" 4
   if [[ ! -f Infos.dat ]]; then
      mv no.Infos.dat Infos.dat
      exit 0
   fi
fi
mkdir calculations
cp $prm.prm calculations/
cp $Project-tk.xyz calculations/$Project.xyz
cp $Project-tk.key calculations/$Project.key
cp $Project-tk.input calculations/$Project.input
cd calculations

# Creation of the directory where the HF optimization will be run and movement of all the files
# If the neutral retinal is used, the right template is copied
#
newdir=${Project}_OptSCF
mkdir ${newdir}
mv ${Project}.key ${newdir}/${newdir}.key
mv ${Project}.xyz ${newdir}/${newdir}.xyz
mv ${prm}.prm ${newdir}/
cp $templatedir/molcas-job.sh ${newdir}/
if [[ $retstereo == "nAT" ]]; then
   cp $templatedir/templateOPTSCFneu ${newdir}/templateOPTSCF 
else
   cp $templatedir/templateOPTSCF ${newdir}/
fi
mv ${Project}.input ${newdir}/${newdir}.input
cd ${newdir}/

# Putting project name, input directory, time and memory in the submission script
# Since it is a SCF optimization, 1500 MB should be enough...
# And 30 hrs are enough as well
#
case $version in
     4.2 | 5.1)
     sed -i "s|MOLCASMEM=MEMORIAMB|MOLCASMEM=MEMORIAMB|" molcas-job.sh
     ;;
     6.2 | 6.3)
     sed -i "s|MOLCASMEM=MEMORIAMB|MOLCAS_MEM=MEMORIAMB|" molcas-job.sh
esac



sed -i "s|NOMEPROGETTO|${newdir}|" molcas-job.sh
no=$PWD
sed -i "s|NOMEDIRETTORI|${no}|" molcas-job.sh
sed -i "s|MEMTOT|1700|" molcas-job.sh
sed -i "s|MEMORIA|1500|" molcas-job.sh
sed -i "s|hh:00:00|30:00:00|" molcas-job.sh

#
# Replacing Set maxiter with MOLCAS_MAXITER for Molcas 8
# Replacing PARAMETRI with current prm filename templateOPTSCF
# 
#
case $version in
     4.2 | 5.1)
     sed -i "s|Set   maxiter   100|Set   maxiter   100|" templateOPTSCF
     ;;
     6.2 | 6.3)
     sed -i "s|Set   maxiter   100|export MOLCAS_MAXITER = 100|" templateOPTSCF
esac

sed -i "s|PARAMETRI|${prm}|" templateOPTSCF

# Generation of the correct Molcas input from templateOPTSCF
#
mv ${newdir}.input temp
cat temp templateOPTSCF > ${newdir}.input
rm temp templateOPTSCF

# Job submission
#
echo " Submitting the HF optimization now..."
echo ""
sleep 1

qsub molcas-job.sh

cd ..
cp $templatedir/Molcami2.sh .
cp $templatedir/finalPDB.sh .

# Updating the Infos.dat Current field, which stores the current running calculation
#
awk '{ if ( $1 == "CurrCalc") sub($2,"OptSCF"); print}' ../Infos.dat > temp
mv temp ../Infos.dat

echo ""
echo "As soon as the HF calculation is finished, cd calculations/ then run Molcami2.sh"
echo ""
