#!/bin/bash
#
# This script performs some analyses on the MD trajectory:
# - plot of RMSD vs time for the atom subset moved during MD
# - plot of RMSD vs time for each residue moved during MD
# - calculation of MD average structure
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" Infos.dat | awk '{ print $2 }'`
backb=`grep "BackBoneMD" Infos.dat | awk '{ print $2 }'`
moving=`grep "MoveSel" Infos.dat | awk '{ print $2 }'`
heatmd=`grep "HeatMD" Infos.dat | awk '{ print $2 }'`
equimd=`grep "EquiMD" Infos.dat | awk '{ print $2 }'`
prodmd=`grep "ProdMD" Infos.dat | awk '{ print $2 }'`
#
# Ask the users about performing the MD analysis
#
echo ""
echo " MD analysis - What do you want to do?"
echo ""
echo " 1) Interactive MD analysis"
echo " 2) Background MD analysis"
echo " 3) No analysis"
echo ""
read canswer
while [[ $canswer -ne 1 &&  $canswer -ne 2 && $canswer -ne 3 ]]; do
      echo " Wrong answer, please type a correct answer"
      read canswer
done
if [[ $canswer -eq 1 || $canswer -eq 2 ]]; then
#
# Check for previous runs of MD analysis and archiving it
#
   if [[ -d Analysis_Dyna ]]; then
      k=1
      while [[ -d $k.Analysis_Dyna ]]; do
            k=$(($k+1))
      done
      mv Analysis_Dyna $k.Analysis_Dyna
   fi
   mkdir Analysis_Dyna
   cd Analysis_Dyna
   cp $templatedir/extraver.sh .
   cp $templatedir/rmsd_residue.sh .
   cp $templatedir/ndx-maker.sh .
   cp ../Dynamic/output/$Project.tpr .
   cp ../Dynamic/output/$Project.gro .
   cp ../Dynamic/output/${Project}.ndx  .
   cp ../Dynamic/output/${Project}.trr .
#
# Generating the selection to be analyzed by calling ndx_maker.sh with the right parameters
# NOTE: the presence of a ndx file in the current folder is assumed
#
   ./ndx-maker.sh $Project 5 $backb y 
#
# Generating the xvg plot with the Gromacs utility g_rms
# Available groups are counted by grepping the square brackets
#
   ngroups=`grep '\[' $Project.ndx | wc -l`
   ngroups=$(($ngroups-1))
   if [[ -f cursel.log ]]; then
      selection=`cat cursel.log`
   else
      echo " No current selection found!"
      echo " Please check what went wrong with ndx-maker.sh"
      echo " Aborting..."
      echo ""
      exit 0
   fi
   echo "$ngroups"> choices.txt
   echo "$ngroups" >> choices.txt
   $gropath/gmx g_rms -s $Project.tpr -f $Project.trr -n $Project.ndx -o rmsd-vs-time.xvg < choices.txt
   rm choices.txt
fi
#
# Plotting rmsd-vs-time.xvg with xmgrace if available, otherwise send message to the user and go on
#
if [[ $canswer -eq 1 ]]; then
   check=`which xmgrace | wc -l`
   if [[ $check -eq 1 ]]; then
      echo ""
      echo " xmgrace found! Press any key to plot the RMSD vs time"
      echo " Then close the window to quit xmgrace and proceed with execution"
      read answer
      echo ""
      xmgrace rmsd-vs-time.xvg
   else
      echo ""
      echo " xmgrace not found! Please plot the rmsd-vs-time.xvg with some other program"
      echo " Going on with execution..."
      echo ""
   fi
fi
#
# Calling the script to plot the RMSD of each residue along the trajectory,
# after computing the beginning and the end of the production phase
#
if [[ $canswer -eq 1 || $canswer -eq 2 ]]; then
   initprod=$(($heatmd+$equimd+1))
   ./rmsd_residue.sh $Project "$selection" $backb $initprod $prodmd
#
# Visualizing all the RMSD plots, after checking for the availability of ImageMagick display command
#
   if [[ $canswer -eq 1 ]]; then
      numimages=`ls -l Images/group*.png | wc -l`
      echo ""
      echo " Showing residues RMSD in groups of four"
      sleep 1
      check=`which display | wc -l`
      if [[ $check -eq 1 ]]; then
         imgview=display
      else
         imgview=gthumb
      fi
      for i in `seq 1 $numimages`; do
          echo " Group $i"
          echo ""
          $imgview Images/group$i.png
      done
   else
      echo ""
      echo " MD analysis completed in background"
      echo " The RMSD per residue plots are stored in Analysis_Dyna/Images"
      echo ""
   fi
#
# Messages to the user
#
   echo " Now Analysis_MD.sh will terminate"
   echo " If MD was OK, run new2mins.sh to go on, otherwise re-run DynIt.sh with a different setup"
   echo ""
else
   echo " You chose not to perform the MD analysis"
   echo " Run new2mins.sh to go on"
   echo ""
fi
exit 0
