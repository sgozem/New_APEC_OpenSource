#!/bin/bash
#
# Reading data from Infos.dat
#
Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
Step=`grep "Step" Infos.dat | awk '{ print $2 }'`
templatedir=`grep "Template" Infos.dat | awk '{ print $2 }'`
chain=`grep "RetChain" Infos.dat | awk '{ print $2 }'`
lysnum=`grep "LysNum" Infos.dat | awk '{ print $2 }'`
prm=`grep "Parameters" Infos.dat | awk '{ print $2 }'`
MultChain=`grep "MultChain" Infos.dat | awk '{ print $2 }'`
solvent=`grep "SolventBox" Infos.dat | awk '{ print $2 }'`
uphess=`grep "Update_hessian" Infos.dat | awk '{ print $2 }'`
moldy=`grep "MD_ensemble" Infos.dat | awk '{ print $2 }'`
updch=`grep "Update_charges" Infos.dat | awk '{ print $2 }'`
ions=`grep "Added_Ions" Infos.dat | awk '{ print $2 }' | awk '{print substr ($0, 0, 1)}'`
amber=`grep "AMBER" Infos.dat | awk '{ print $2 }'`

if [ -d ../Step_$(($Step+1)) ]; then
      echo " Folder \"Step_$(($Step+1))\" found! Something is wrong ..."
      echo " Terminating ..."
      exit 0
      echo ""
fi

mkdir ../Step_$(($Step+1))
mkdir ../Step_$(($Step+1))/Dynamic
cp Infos.dat ../Step_$(($Step+1))
./update_infos.sh "Step" $(($Step+1)) ../Step_$(($Step+1))/Infos.dat
cp update_infos.sh $prm.prm template_* ../Step_$(($Step+1))
cp -r Chromophore ../Step_$(($Step+1))
cp $templatedir/ASEC/MD_NVT.sh ../Step_$(($Step+1))
cp $templatedir/ASEC/dynamic_sol_NVT.mdp ../Step_$(($Step+1))/Dynamic
cp $templatedir/gromacs.sh ../Step_$(($Step+1))/Dynamic
cp -r $templatedir/$amber.ff ../Step_$(($Step+1))/Dynamic
#cp ../Step_$(($Step+1))/Dynamic/$amber.ff/normalamino-h ../Step_$(($Step+1))/Dynamic/$amber.ff/aminoacids.hdb
#cp ../Step_$(($Step+1))/Dynamic/$amber.ff/amino-rettrans ../Step_$(($Step+1))/Dynamic/$amber.ff/aminoacids.rtp
cat calculations/ESPF_charges/new_rtp >> ../Step_$(($Step+1))/Dynamic/$amber.ff/aminoacids.rtp  
cp calculations/new_charges ../Step_$(($Step+1))
cp calculations/${Project}_VDZP_Opt/${Project}_VDZP_Opt.JobIph_new ../Step_$(($Step+1))/${Project}_VDZP_Opt.JobIph_old
cp calculations/${Project}_finalPDB/${Project}_new.gro ../Step_$(($Step+1))/Dynamic/${Project}_box_sol.gro
cp Dynamic/${Project}_box_sol.top ../Step_$(($Step+1))/Dynamic
cp Dynamic/${Project}_box_sol.ndx ../Step_$(($Step+1))/Dynamic
cp Dynamic/residuetypes.dat ../Step_$(($Step+1))/Dynamic
cp Dynamic/*.itp ../Step_$(($Step+1))/Dynamic

echo ""
echo " *******************************************************"
echo ""    
echo " Go to Step_$(($Step+1)) and continue with \"MD_NVT.sh\""
echo ""
echo " *******************************************************"

