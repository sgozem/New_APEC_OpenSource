#!/bin/bash
#
# This installation script is supposed to reside in the same folder as template and New_APEC.sh
# resulting from the github download. 
#
cd template

# Now the user interface starts
#
echo " Welcome to the installation of APEC (Average Protein Electrostatic Configuration)"
echo ""
echo " brought to you by Gozem Lab"
echo " Georgia State University"
echo " Developed by Yoelvis Orozco-Gonzalez"
echo " and contributions from Federico Melaccio"
echo ""

# It is assumed that the current template dir will be permanent,
# therefore the Infos.dat is created with that information if user confirms
#
echo " The current directory $PWD will be the source of the template files"
echo " Please type y to confirm or n to stop installation"
read answer
echo ""
while [[ $answer != y && $answer != n ]]; do
	echo " Wrong choice, please type y or n"
	read answer
	echo ""
done
if [[ $answer == n ]]; then
   echo " Installation is aborting. Goodbye!"
   echo ""
   exit 0
fi

echo "Template $PWD" > Infos.dat
templatepath=$PWD
echo ""

# Checking if a configuration file is available
# If so, all the information is read from there, otherwise
# the user is asked for the missing data
#
if [[ -f ../APECconfig ]]; then
   echo " APECconfig configuration file was found! Do you want to use it? (y/n)"
   read choiceconfig
   echo ""
   while [[ $choiceconfig != y && $choiceconfig != n ]]; do
         echo " Wrong choice! Please type either y or n!"
         read choiceconfig
         echo ""
   done
else
   choiceconfig=n
fi

# Asking for Linux commands absolute path
# This is required for portability
#
echo " Now you have to provide some data to implement APEC"
echo ""
echo " 1) Linux Commands"
echo ""
echo "    APEC needs to know which folder includes the Linux commands"
echo ""
echo "    NOTE: in some Unix or Linux implementation (like AIX),"
echo "          a special folder includes the Linux-like commands,"
echo "          e.g. /usr/linux/bin"
echo ""
echo "    If you don't have any clue, type"
echo ""
echo "    which ls"
echo ""
echo "    and see what is the path preceeding ls, e.g. if you see /usr/bin/ls"
echo "    type here /usr/bin (last slash must be omitted)."
echo "    If you are not sure, type NO"
read pathway
if [[ $pathway == NO ]]; then
   echo ""
   echo "    You decided to use the standard Linux commands"
   echo ""
   sed -i "s|COMMANDER/||g" *.sh
   linuxpath=""
   pathway=""
   commander="COMMANDER/"
else
   echo ""
   echo "    You typed $pathway, updating all the APEC files..."
   echo ""
   linuxpath=$pathway/
   "$linuxpath"sed -i "s|COMMANDER|$pathway|g" *.sh
   commander="COMMANDER"
fi

# Changing INSTALLDIR in New_APEC.sh to have the path for Infos.dat
#
cd ..
"$linuxpath"sed -i "s|INSTALLDIR|$PWD|g" New_APEC.sh
"$linuxpath"sed -i "s|$commander|$pathway|g" New_APEC.sh
cd template

# Now the absolute path of Molcas installation
#
if [[ $choiceconfig == y ]]; then
   pathway=`"$linuxpath"grep 'MolcasInstall' ../APECconfig | awk '{ print $2 }'`
else
   echo " 2) Molcas installation directory"
   echo ""
   echo "    Please type the absolute path of the Molcas installation you wish to use"
   echo "    Omit the last slash of the path"
   read pathway
   echo ""
fi
if [[ -f $pathway/Symbols ]]; then
   echo "    You will use Molcas in $pathway, updating all the APEC files..."
   echo ""
   "$linuxpath"sed -i "s|MOLCASDIR|$pathway|g" molcas*
else
   echo "    Wrong molcas folder, installation is aborting..."
   echo ""
   rm Infos.dat
   exit 0
fi

# Now the absolute path of Molcas driver is requested
#
if [[ $choiceconfig == y ]]; then
   pathway=`"$linuxpath"grep 'MolcasDriver' ../APECconfig | awk '{ print $2 }'`
else
   echo " 3) Molcas driver"
   echo ""
   echo "    Please type the absolute path of Molcas execution script"
   echo "    Usually it is stored in a bin/ directory somewhere"
   echo "    Omit the last slash of the path"
   read pathway
   echo ""
fi
if [[ -f $pathway/molcas ]]; then
   echo "    You will use Molcas in $pathway, updating all the APEC files..."
   echo ""
   "$linuxpath"sed -i "s|MOLCASDRV|$pathway|g" molcas*
else
   echo "    molcas driver not found, installation is aborting..."
   echo ""
   rm Infos.dat 
   exit 0
fi

# Asking the Molcas version to select the Espf input format
#
if [[ $choiceconfig == y ]]; then
   version=`"$linuxpath"grep 'MolcasVersion' ../APECconfig | awk '{ print $2 }'`
   case $version in
        7.2 | 7.4 | 7.5) truever=1
        ;;
        7.6 | 7.7 | 7.8 | 7.9) truever=2
        ;;
        8.0 | 8.1 ) truever=3               
        ;;                                  
        *) truever=4                        
   esac
   version=$truever
else
   echo " 4) Molcas version"
   echo ""
   echo "    Please select the Molcas version you will be using"
   echo "    1. - version 7.5 or below"
   echo "    2. - version 7.6 - 7.9"         
   echo "    3. - version 8.0 or above"      
   read version
fi
case $version in
	1) 
	"$linuxpath"sed -i "s|External = Tinker|External = @tinker|g" template*
	echo "    ESPF input will include External = @tinker"
        echo ""
	;;
	2)
	echo "    ESPF input will include External = Tinker"
	echo ""
	;;
	3)
	echo "    ESPF input will include External = Tinker"
	echo ""
	;;
	*)
	echo "    Wrong selection, installation is aborting..."
	echo ""
	rm Infos.dat
	exit 0
esac

# Now the absolute path of Tinker executable is requested
#
if [[ $choiceconfig == y ]]; then
   pathway=`"$linuxpath"grep 'TinkerDir' ../ARMconfig | awk '{ print $2 }'`
else
   echo " 5) Tinker executables directory"
   echo ""
   echo "    Please type the absolute path where the Tinker executables files are stored"
   echo "    It is usually the bin/ folder in the Tinker installation directory"
   echo "    Omit the last slash of the path"
   read pathway
   echo ""
fi
if [[ -f $pathway/pdbxyz ]]; then
   echo "    You will use the Tinker exe files in $pathway, updating the ARM files..."
   echo ""
   echo "Tinker $pathway" >> Infos.dat
   "$linuxpath"sed -i "s|TINKERDIR|$pathway|g" molcas*
else
   echo "    Tinker executables not found, installation is aborting..."
   echo ""
   rm Infos.dat 
   exit 0
fi

# Asking the Tinker version to know the right file formats
#
if [[ $choiceconfig == y ]]; then
   version=`"$linuxpath"grep 'TinkerVersion' ../ARMconfig | awk '{ print $2 }'`
   case $version in
        4.2) truever=1
        ;;
        5.1) truever=2
        ;;
        6.2 | 6.3 ) truever=3
        ;;
        *) truever=4
   esac
   version=$truever
else
   echo " 6) Tinker version"
   echo ""
   echo "    Please select the Tinker version you will be using"
   echo "    1. - version 4.2"
   echo "    2. - version 5.1 - 5.9"
   echo "    3. - version 6.3 or above"
   read version
fi
case $version in
        1)
        echo "Version Tk 4.2" >> Infos.dat
        ;;
        2)
        echo "Version Tk 5.1" >> Infos.dat
        ;;
        3)
        echo "Version Tk 6.3" >> Infos.dat
        ;;
        4)
        echo " Wrong selection, installation is aborting..."
        echo ""
        rm Infos.dat
        exit 0
esac


# Patching the original Tinker 4.2 or 5.1 or 6.3 code for retinal recognition
# and other options like adding link atoms
# If the utility patch is not existing, the whole section is skipped
#
rigapatch=`which "$linuxpath"patch | "$linuxpath"head -n 1`
testpatch=${rigapatch:0:1}
if [ $testpatch == "/" ]; then
   cd $pathway
   cd ..
   if [ -d source/ ]; then
   	cd source/
        checkret=`grep 'buildret' pdbxyz.f | head -n 1`
        if [[ -z $checkret ]]; then
           cd ..
   	   cp -r source/ source-orig/ 
   	   echo "    The source/ directory in $PWD was backed up as source-orig/"
   	   echo "    Preparing to patch Tinker source code..."
   	   echo ""
           cd source/
           case $version in
                1)
                "$linuxpath"cp $templatepath/tk-4.2-rhod.diff .
                "$linuxpath"patch < tk-4.2-rhod.diff
                ;;
                2)
                "$linuxpath"cp $templatepath/rhodopsin-tk.diff .
                "$linuxpath"patch < rhodopsin-tk.diff
                ;;
                3)
                "$linuxpath"cp $templatepath/rhodopsin-tk6.diff .
                "$linuxpath"patch < rhodopsin-tk6.diff
                ;;
                *)
                echo " Wrong selection, installation is aborting..."
                echo ""
                rm Infos.dat
                exit 0
           esac
	   checkaci=`grep 'buildaci' pdbxyz.f | head -n 1`
           if [[ -z $checkaci ]]; then
              echo "    Adding support for free acetate ion (ACI)..."
              echo ""
              "$linuxpath"cp $templatepath/pdbxyz.diff .
              "$linuxpath"patch < pdbxyz.diff
           else
              echo "    WARNING: Free acetate ion code detected!"
              echo "    ARM will not modify Tinker code for that"
              echo ""
           fi
   	   echo ""
           echo "    Tinker source code patched successfully!"
   	   echo "    Now recompiling it..."
   	   echo ""
   	   "$linuxpath"make all > make.log
   	   error=`grep "Error" make.log | wc -l`
           if [ $error -ne 0 ]; then
                echo "    Tinker compilation went wrong! Installation is aborting..."
   	        echo ""
   	        cd $templatepath/template
   	        rm Infos.dat
   	        exit 0
           else
   	        rm make.log
                "$linuxpath"make rename
   	        echo ""
   	        echo "    Tinker compiled successfully!"
   	        echo ""
           fi
        else
           echo "    WARNING: Retinal building code detected!"
           echo "    ARM will assume that Tinker source code is already patched."
           echo ""
           checkaci=`grep 'buildaci' pdbxyz.f | head -n 1`
           if [[ -z $checkaci ]]; then
              echo "    Adding support for free acetate ion (ACI)..."
              echo ""
              "$linuxpath"cp $templatepath/pdbxyz.diff .
              "$linuxpath"patch < pdbxyz.diff
              "$linuxpath"make all > make.log
              error=`grep "Error" make.log | wc -l`
              if [ $error -ne 0 ]; then
                 echo "    Tinker compilation for free acetate went wrong! Installation is aborting..."
                 echo ""
                 cd $templatepath/template
                 rm Infos.dat
                 exit 0
              else
                 rm make.log
                 "$linuxpath"make rename
                 echo ""
                 echo "    ACI support added successfully!"
                 echo ""
              fi
           else
              echo "    WARNING: Free acetate ion code detected!"
              echo "    ARM will not modify Tinker code for that"
              echo ""
           fi
        fi
   else
   	echo "    WARNING: Tinker source directory not found!"
   	echo "    You will need to manually patch Tinker source code,"
   	echo "    otherwise retinal will not be recognized..."
   	echo ""
   fi
else
   echo "   WARNING: patch utility not found in $linuxpath!"
   echo "   You will need to manually patch Tinker source code,"
   echo "   otherwise retinal will not be recognized..."
   echo ""
fi
cd $templatepath

# Now the type of queuing system and the associated command
#
if [[ $choiceconfig == y ]]; then
   answer=`"$linuxpath"grep 'Queue' ../ARMconfig | awk '{ print $2 }'`
   case $answer in
        LoadLeveler) realans=1
        ;;
        Torque) realans=2
        ;;
        SunGridEngine) realans=3
        ;;
	OAR) realans=4
	;;
	SLURM) realans=5
	;;
        *) realans=6
        ;;
   esac
   answer=$realans
else
   echo " 7) Queue system and command"
   echo ""
   echo "    According to the queue system of the cluster you are using,"
   echo "    the submission command and scripts change accordingly"
   echo "    Please choose the appropriate option:"
   echo "    1. - LoadLeveler, on AIX (molcas-tk.cmd, command: llsubmit)"
   echo "    2. - Torque/PBS, on OSC Glenn cluster (molcas-job.sh, command: qsub)"
   echo "    3. - SunGridEngine, on all Rocks cluster (molcas.sub, command: qsub)"
   echo "    4. - OAR, on Rheticus/Marseille cluster (molcas.oar.sh, command: oarsub -S)"
   echo "    5. - SLURM, GSU and comet clusters (molcas.icpms.sh, command: sbatch)"
   echo "    6. - None of the above"
   read answer
   echo ""
   while [[ $answer -gt 6 || $answer -lt 1 ]]; do
	 echo " Wrong choice, please type a number between 1 and 5"
	 read answer
	 echo ""
   done
fi
case $answer in
	1)
	"$linuxpath"sed -i "s|SUBMISSION|molcas-tk.cmd|g" *.sh
	"$linuxpath"sed -i "s|MINSUB|minimize.cmd|g" *.sh
	"$linuxpath"sed -i "s|DYNSUB|gromacs.cmd|g" *.sh
	"$linuxpath"sed -i "s|SUBCOMMAND|llsubmit|g" *.sh
#	"$linuxpath"sed -i "s|QUEUECHECK|llq|g" *
#	"$linuxpath"sed -i "s|QCHAR|I|g" *
#	"$linuxpath"sed -i "s|MODULELOAD|module load|g" *
	echo "Queue AIX" >> Infos.dat
	;;
	2)
	"$linuxpath"sed -i "s|SUBMISSION|molcas-job.sh|g" *.sh
	"$linuxpath"sed -i "s|MINSUB|minimize.sh|g" *.sh
	"$linuxpath"sed -i "s|DYNSUB|gromacs.sh|g" *.sh
        "$linuxpath"sed -i "s|SUBCOMMAND|qsub|g" *.sh
#	"$linuxpath"sed -i "s|QUEUECHECK|qstat|g" *
#	"$linuxpath"sed -i "s|QCHAR|Q|g" *
#	"$linuxpath"sed -i "s|MODULELOAD vmd||g" *
	echo "Queue Glenn" >> Infos.dat
	;;
	3)
	"$linuxpath"sed -i "s|SUBMISSION|molcas.sub|g" *.sh
        "$linuxpath"sed -i "s|MINSUB|minimize.sub|g" *.sh
        "$linuxpath"sed -i "s|DYNSUB|gromacs.sub|g" *.sh
	"$linuxpath"sed -i "s|SUBCOMMAND|qsub|g" *.sh
#	"$linuxpath"sed -i "s|QUEUECHECK|qstat|g" *
#       "$linuxpath"sed -i "s|QCHAR|qw|g" *
#	"$linuxpath"sed -i "s|MODULELOAD vmd||g" *
        echo "Queue Rocks" >> Infos.dat	
	;;
	4)
	"$linuxpath"sed -i "s|SUBMISSION|./molcas.oar.sh|g" *.sh
	"$linuxpath"sed -i "s|MINSUB|./minimize.oar.sh|g" *.sh
	"$linuxpath"sed -i "s|DYNSUB|./gromacs.oar.sh|g" *.sh
	"$linuxpath"sed -i "s|SUBCOMMAND|oarsub -S|g" *.sh
	echo "Queue Rheticus" >> Infos.dat 
	;;
	5)
        "$linuxpath"sed -i "s|SUBMISSION|molcas.slurm.sh|g" *.sh
        "$linuxpath"sed -i "s|SUBMIGAU|gaussian.slurm.sh|g" *.sh
        "$linuxpath"sed -i "s|MINSUB|minimize.slurm.sh|g" *.sh
        "$linuxpath"sed -i "s|DYNSUB|gromacs.slurm.sh|g" *.sh
        "$linuxpath"sed -i "s|SUBCOMMAND|sbatch|g" *.sh
        "$linuxpath"sed -i "s|SUBMISSION|molcas.slurm.sh|g" ASEC/*.sh
        "$linuxpath"sed -i "s|SUBMIGAU|gaussian.slurm.sh|g" ASEC/*.sh
        "$linuxpath"sed -i "s|MINSUB|minimize.slurm.sh|g" ASEC/*.sh
        "$linuxpath"sed -i "s|DYNSUB|gromacs.slurm.sh|g" ASEC/*.sh
        "$linuxpath"sed -i "s|SUBCOMMAND|sbatch|g" ASEC/*.sh

        echo "Queue SLURM" >> Infos.dat	
	;;
	6)
	echo " Your cluster and queue system is not implemented yet"
	echo " Installation is aborting..."
	echo ""
	rm Infos.dat
	exit 0
	;;
	*)
	echo " Something went wrong. Installation is aborting..."
	echo ""
	rm Infos.dat
	exit 0
esac

# Asking if Gromacs will be needed, and the path of its executables
#
if [[ $choiceconfig == y ]]; then
   checkgro=`"$linuxpath"grep 'GromacsPath' ../ARMconfig | awk '{ print $2 }'`
   if [[ -z $checkgro ]]; then
      gromacs=n
   else
      gromacs=y
   fi
else
   echo " 8) Gromacs - testing"
   echo ""
   while [[ $gromacs != "y" && $gromacs != "n" ]]; do
         echo "    Please type y/n if you want to use Gromacs"
         read gromacs
   done
fi
case $gromacs in
        y)
        echo "Gromacs YES" >> Infos.dat
        if [[ $choiceconfig == y ]]; then
           pathgro=$checkgro
        else
           echo "    Please type the Gromacs executables path you wish to use"
           read pathgro
        fi
        echo "GroPath $pathgro" >> Infos.dat
        ;;
        n)
        echo "Gromacs NO" >> Infos.dat
        ;;
        *)
        echo " Wrong selection, installation is aborting..."
        echo ""
        rm Infos.dat
        exit 0
esac

# Asking if DOWSER will be needed, and the path of its executables
# DOWSER is a program for water and polar hydrogen placement
#
if [[ $choiceconfig == y ]]; then
   dowserchk=`"$linuxpath"grep 'DowserPath' ../ARMconfig | awk '{ print $2 }'`
   if [[ -z $dowserchk ]]; then
      dowser=n
   else
      dowser=y
   fi
else
   echo " 9) Dowser"
   echo ""
   while [[ $dowser != "y" && $dowser != "n" ]]; do
         echo "    Please type y/n if you want to use Dowser"
         read dowser
   done
fi
case $dowser in
        y)
        echo "Dowser YES" >> Infos.dat
        if [[ $choiceconfig == y ]]; then
           pathdow=$dowserchk
        else
           checkdow=`which dowser | wc -l`
           pathdow=`which dowser`
           while [[ $checkdow -ne 1 ]]; do
                 echo "    dowser executable not found in the system path!"
                 echo "    Please type the full Dowser path you wish to use"
                 read pathdow
                 if [[ -f $pathdow ]]; then
                    checkdow=1
                 else
                    checkdow=0
                 fi  
           done
        fi
        echo "DowPath $pathdow" >> Infos.dat
        ;;
        n)
        echo "Dowser NO" >> Infos.dat
        ;;
        *)
        echo " Wrong selection, installation is aborting..."
        echo ""
        rm Infos.dat
        exit 0
esac

# Asking for mutants stuff 
#
#
if [[ $choiceconfig == y ]]; then
   modelconf=`"$linuxpath"grep 'Modeller' ../ARMconfig | awk '{ print $2 }'`
   scwrlconf=`"$linuxpath"grep 'SCWRL4' ../ARMconfig | awk '{ print $2 }'`
   if [[ $modelconf == YES ]]; then
      chkmodel=`which mod9.11`
      if [[ -z $chkmodel ]]; then
         chkmodel=`which mod9.12`
      fi
      if [[ -z $chkmodel ]]; then
         modelconf=NO
      fi
   fi
   echo "Modeller $modelconf" >> Infos.dat
   if [[ $scwrlconf == YES ]]; then
      chkscwrl=`which Scwrl4`
      if [[ -z $chkscwrl ]]; then
         scwrlconf=NO
      fi
   fi
   echo "SCWRL4 $modelconf" >> Infos.dat
else
   echo " 10) Mutations - Modeller and SCWRL4"
   echo ""
   mutations="s"
   while [[ $mutations != "y" && $mutations != "n" ]]; do
         echo "    Do you want to use Modeller and/or SCWRL4 for mutations? (y/n)"
         read mutations
   done
   if [[ $mutations == n ]]; then
      echo "    No mutants will be possible with the current installation"
      echo ""
      echo "Modeller NO" >> Infos.dat
      echo "SCWRL4 NO" >> Infos.dat
   else
      chkmodel=`which mod9.11`
      if [[ -z $chkmodel ]]; then
         chkmodel=`which mod9.12`
      fi
      chkscwrl=`which Scwrl4`
      echo ""
      if [[ -z $chkmodel && -z $chkscwrl ]]; then
         echo " Neither Modeller nor SCWRL4 were found in the system path!"
         echo " Please be sure you installed them properly"
         echo " No mutants will be possible with the current installation"
         echo ""
         echo "Modeller NO" >> Infos.dat
         echo "SCWRL4 NO" >> Infos.dat
      else
         if [[ -z $chkmodel ]]; then
            echo " Modeller was not found in the system path!"
            echo " If you installed it, please be sure it is in your PATH and restart the installation"
            echo " Just SCWRL4 will be used, as it is in $chkscwrl"
            echo ""
            echo "Modeller NO" >> Infos.dat
            echo "SCWRL4 YES" >> Infos.dat
         else
            if [[ -z $chkscwrl ]]; then
               echo " SCWRL4 was not found in the system path!"
               echo " If you installed it, please be sure it is in your PATH and restart the installation"
               echo " Just Modeller will be used, as it is in $chkmodel"
               echo ""
               echo "Modeller YES" >> Infos.dat
               echo "SCWRL4 NO" >> Infos.dat
            else
               echo " Modeller was found in $modeller"
               echo ""
               echo "Modeller YES" >> Infos.dat
               echo " SCWRL4 was found in $scwrl"
               echo ""
               echo "SCWRL4 YES" >> Infos.dat
            fi
         fi
      fi
   fi
fi
#
# When Modeller is found, the user has to check if the current system needs to be told about Modeller
# in LD_LIBRARY_PATH and PYTHONPATH. Otherwise it won't work
#
if [[ -n $chkmodel ]]; then
   echo " NOTE: if Modeller is not working, you might need to update LD_LIBRARY_PATH and PYTHONPATH as follows:"
   echo ' export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:[path-to-modlib]:[path-to-lib/arch]'
   echo ""
   echo ' export PYTHONPATH=$PYTHONPATH:[path-to-modlib]:[path-to-lib/arch]'
   echo ""
   echo " where the two folders are in the Modeller installation folder"
   echo ""
fi

# Successful installation message to the user 
#
echo " All the required information was gathered"
echo " ARM was installed successfully"
echo ""
echo " Please refer to documentation to start working"
echo " Have a nice day!"
echo ""
