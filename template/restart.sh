#!/bin/bash
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../../Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" ../../Infos.dat | awk '{ print $2 }'`
tinkerdir=`grep "Tinker" ../../Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" ../../Infos.dat | awk '{ print $2 }'`
version=`grep "Tk" ../../Infos.dat | awk '{ print $3 }'`
# Instructions to the user
#
echo ""
echo " The current project is $Project. Which optimization do you want to continue?"
echo " Type 1 for 3-21G and 2 for 6-31G*"
read answer
echo ""
while [[ $answer != 1 && $answer != 2 ]]; do
      echo ""
      echo " Wrong selection. Please type 1 for 3-21G or 2 for 6-31G*"
      read answer
      echo ""
done

# Deciding on the name of the folder
#
if  [ $answer == 1 ]; then
	folder="$Project"_3-21G_Opt
else
	folder="$Project"_6-31G_Opt
fi

# If the directory is not found, the script is terminated
#
if [ ! -d $folder ]; then
	echo " Directory $folder not found. Please check if everything is ok."
	echo " restart.sh is terminating..."
	echo ""
	exit 0
fi

# Searching for previous run of the same calculation and renaming folders
#
number=`ls -ld * | grep "$folder" | wc -l | awk '{ print $1 }'`
if [ $number == 0 ]; then
	echo " This is the first time you restart this calculation"
	echo " Files from the previous run will be moved to $folder.1"
	echo ""
	mv $folder $folder.1
	lastdir=$folder.1
else
	echo " This optimization was run $number times"
	echo ""
	i=1
	dir=$folder.$i
	while [ -d $dir ]; do
		i=$(($i+1))
		dir=$folder.$i
	done
	echo " Files of the last run are moved in $folder.$i"
	echo ""
        mv $folder $folder.$i
	lastdir=$folder.$i
fi

# Preparing the files and folder for the new calculation. 
#
mkdir $folder
cp $lastdir/$folder.Final.xyz $folder/$folder.xyz
cp $lastdir/$folder.key $folder/
cp $lastdir/$prm.prm $folder/
cp $lastdir/molcas-job.sh $folder/
cp $lastdir/$folder.JobIph $folder/
cp $lastdir/$folder.Espf.Data $folder/

# modify-inp.vim and the template for CASSCF optimization are copied from $templatedir
#
cp $templatedir/modify-inp.vim $folder/
cp $templatedir/template $folder/

# The input file must be generated again with xyzedit on Final.xyz from last run
#
cd $folder
case $version in
     4.2 | 5.1)
     echo "20" >> tinkerxyzedit.vim
     echo "1" >> tinkerxyzedit.vim
     ;;
     6.2 | 6.3)
     echo "21" >> tinkerxyzedit.vim
     echo "1" >> tinkerxyzedit.vim
esac

$tinkerdir/xyzedit $folder.xyz < tinkerxyzedit.vim
rm tinkerxyzedit.vim


# Editing the input file with modify-inp.vim: removal of the Tinker standard part
# and introduction of the basis set, which depends on $answer variable
#
if [ $answer == 1 ]; then
	sed -i 's/BASIS/3-21G/' modify-inp.vim
        sed -i 's/BAS2/3-21G/' modify-inp.vim
else
	sed -i 's/BASIS/6-31G*/' modify-inp.vim
	sed -i 's/BAS2/6-31G\\*/' modify-inp.vim
fi
vim -es $folder.input < modify-inp.vim
rm modify-inp.vim

# Modifying the template for CAS opt to use the correct parameter file
#
sed -i "s|PARAMETRI|${prm}|" template

# Merging the geometrical part and the template for CAS optimization
#
mv $folder.input temp
cat temp template > $folder.input
rm temp template

# The calculation is ready to be submitted
#
echo ""
echo " Restarting the $folder calculation now..."
echo ""
sleep 1

qsub molcas-job.sh

echo ""
echo " After the completion of this restarted calculation,"
echo " you can go on with the scripts flow as usual"
echo " restart.sh is terminating normally"
echo ""

