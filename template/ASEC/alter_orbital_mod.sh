#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../Infos.dat | awk '{ print $2 }'`
currcalc=`grep "CurrCalc" ../Infos.dat | awk '{ print $2 }'`

# Instructions to the user
#
echo ""
echo " Your project is $Project. The last calculation is $currcalc"
echo " Do you want to change the orbitals from $currcalc? (y/n)"
read answer
echo ""

# If the answer is not yes, the script is terminated
#
if [ $answer == "y" ]; then
	echo " Proceeding with orbital switch..."
	echo ""
else
	echo " No orbital switch will be perfomed"
	echo " alter_orbital.sh is terminating..."
	echo ""
	exit 0
fi

# Orbitals are changed according to the current calculation
# The folder is renamed with a tag and the orbital switch is performed
#
case $currcalc in
	VDZ)
#	echo " Orbitals from CAS/VDZ single point will be changed"
#	echo ""
	folder="$Project"_VDZ
	;;
	VDZP)
#	echo " Orbitals from CAS/VDZP single point will be changed"
#	echo ""
	folder="$Project"_VDZP
	;;
	*)
	echo " Orbitals cannot be changed at this step!"
	echo " If they are wrong after an optimization, you might need to re-run it"
	echo " alter_orbital.sh is terminating..."
	exit 0
esac

# The folder is renamed as the first available one with the "wrong" tag
# to allow the execution of orbital switching multiple times
# The orbital switching is always executed in $folder
#
i=1
while [ -d $folder.wrong.$i ]; do
	i=$(($i+1))
done
mv $folder $folder.wrong.$i
mkdir  $folder

cd $folder.wrong.$i
cp $folder.input $folder.key *.xyz *.prm $folder.RasOrb *.sh ../$folder
cd ../$folder

sed -i "s,*>>> COPY \$InpDir/\$Project.RasOrb INPORB,>>> COPY \$InpDir/\$Project.RasOrb INPORB," $folder.input
sed -i "s/Inactive=/* Inactive=/" $folder.input
sed -i "s/Ras2=/* Ras2=/" $folder.input


# Job submission and template copy for the following step
#
echo ""
echo ""
echo " Please go to the folder \"$folder\" and modify the bottom of the"
echo " file \"$Project.RasOrb\" to redifine the new active space"
echo " and resubmit again using \"sbatch molcas-job.sh\" "
echo ""
echo ""
#
case $currcalc in
	VDZ)
	echo " As soon as the calculation is finished, run 1st_to_2nd_mod.sh"
	echo ""
	;;
	VDZP)
	echo " As soon as the calculation is finished, run sp_to_opt_VDZP.sh"
	echo ""
	;;
esac

