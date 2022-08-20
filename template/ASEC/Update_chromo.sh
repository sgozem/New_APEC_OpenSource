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

echo ""
echo " Please wait ..."
echo ""

if [[ $solvent == "YES" ]]; then
   numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
   protatoms=`grep "protatoms" Infos.dat | awk '{ print $2 }'`
else
   numatoms=`grep "numatoms" Infos.dat | awk '{ print $2 }'`
   protatoms=$numatoms
fi

if [[ $uphess ==  "YES" ]]; then
   fromhess=`grep "Hessian_from" Infos.dat | awk '{ print $2 }'`
fi

folder='6-31G_Opt'

mkdir Update_chromo
cd Update_chromo

cp ../calculations/${Project}_finalPDB/${Project}_final.xyz final.xyz
cp ../calculations/${Project}_$folder/$prm.prm .
cp ../template_tk2gro .
cp ../template_gro2tk .
if [[ $solvent == "YES" ]]; then
   cp ../MD_ASEC/Best_Config_full.gro Best_Config.gro
else
   if [[ $Step -eq 0 ]]; then
      cp ../Dynamic/$Project.gro Best_Config.gro
   else
      cp ../MD_ASEC/Best_Config.gro .
   fi
fi

if [[ $updch ==  YES ]]; then
     
   cp ../calculations/RESP_charges/new_charges .
   cp $templatedir/ASEC/Update_chromo.f .

   if [[ $solvent ==  "YES" ]]; then
      Project="${Project}_box_sol"
   fi
#   if [[ $MultChain ==  "NO" ]]; then
#         cp ../Dynamic/${Project}.top chain.itp_old
#   else
#      cp ../Dynamic/${Project}_Protein_chain_$chain.itp chain.itp_old
#   fi
   if [[ -f ../Dynamic/${Project}_Protein_chain_$chain.itp ]]; then
      cp ../Dynamic/${Project}_Protein_chain_$chain.itp chain.itp_old
   else
      cp ../Dynamic/${Project}.top chain.itp_old
   fi
   if [[ $solvent ==  "YES" ]]; then
      Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
   fi
else
   cp $templatedir/ASEC/Update_chromo.f .
   sed -i "s/c\ \ \ \ \ \ GO TO 20/\ \ \ \ \ \ GO TO 20/" Update_chromo.f
   sed -i "s/      open(7/c      open(7/" Update_chromo.f
   sed -i "s/      open(8/c      open(8/" Update_chromo.f
   sed -i "s/      open(9/c      open(9/" Update_chromo.f
fi

##########
# Reading the number of QM and MM atoms from .key file and copy them to the Update_chromo.f
# and to the Charges.f
##########

sed -i "s|numero|$numatoms|g" Update_chromo.f
if [[ $solvent == "YES" ]]; then
   sed -i "s|proatoms|$(($protatoms+$ions))|g" Update_chromo.f
else
   sed -i "s|proatoms|$numatoms|g" Update_chromo.f
fi

interv_mm=`grep -w -c "MM" ../calculations/${Project}_$folder/${Project}_$folder.key`
interv_qm=`grep -w -c "QM" ../calculations/${Project}_$folder/${Project}_$folder.key`
intervals=$(($interv_mm+$interv_qm))

sed -i "s|intervalos_qm|$interv_qm|g" Update_chromo.f
sed -i "s|intervalos|$intervals|g" Update_chromo.f


for i in $(eval echo "{$interv_qm..1}")
do 
   final=`grep -w -m $i "QM" ../calculations/${Project}_$folder/${Project}_$folder.key | tail -n1 | awk '{print $3}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ iqmintv($i,2)=$final" Update_chromo.f
   init=`grep -w -m $i "QM" ../calculations/${Project}_$folder/${Project}_$folder.key | tail -n1 | awk '{print ((-1*$2))}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ iqmintv($i,1)=$init" Update_chromo.f
done

for i in $(eval echo "{$intervals..1}")
do 
   temp=`grep -w -m $i "MM\|QM" ../calculations/${Project}_$folder/${Project}_$folder.key | tail -n1 | awk '{print $3}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,2)=$temp" Update_chromo.f
   temp=`grep -w -m $i "MM\|QM" ../calculations/${Project}_$folder/${Project}_$folder.key | tail -n1 | awk '{print ((-1*$2))}'`
   sed -i "/CCCCCCCCC  Data/a\ \ \ \ \ \ ichromo($i,1)=$temp" Update_chromo.f
done

if [[ $updch == YES ]]; then

   temp1=`grep -A1 "residue $lysnum" chain.itp_old | tail -n1 | awk '{print $3}'`
   temp2=`grep -A1 "residue $lysnum" chain.itp_old | tail -n1 | awk '{print $4}'`

   numfix=`grep -c "$temp1    $temp2" chain.itp_old`
   lastg=`grep -w -m $numfix "${temp1}${temp2}" Best_Config.gro | tail -n1 | awk '{print $3}'`
   sed -i "s|ultimo|$lastg|g" Update_chromo.f
   firstg=`grep -w -m 1 "${temp1}${temp2}" Best_Config.gro | tail -n1 | awk '{print $3}'`
   sed -i "s|primero|$firstg|g" Update_chromo.f
   linefirst=`grep -n "; residue ${temp1} ${temp2}" chain.itp_old | awk '{print $1}' FS=":"`

   fileend=`wc -l chain.itp_old | awk '{print $1}'`

   head -n $linefirst chain.itp_old > part1
   head -n $(($linefirst+$lastg-$firstg+1)) chain.itp_old | tail -n $(($lastg-$firstg+1)) > part2_old
   carbcd=`grep "RET     CD" part2_old | awk '{ print $1 }'`
   sed -i "s/CDNUMERO/$carbcd/" Update_chromo.f
   tail -n $(($fileend-($linefirst+$lastg-$firstg+1))) chain.itp_old > part3

fi

gfortran Update_chromo.f -o Update_chromo.x
./Update_chromo.x

if [[ $updch == YES ]]; then
   cat part1 part2_new part3 > a
   mv a chain.itp
fi

if [[ $solvent ==  "YES" ]]; then
   Project="${Project}_box_sol"
fi

mkdir ../../Step_$(($Step+1))
mkdir ../../Step_$(($Step+1))/Dynamic
cp ../Infos.dat ../../Step_$(($Step+1))
cp ../update_infos.sh ../../Step_$(($Step+1))
sed -i "s|Step $Step|Step $(($Step+1))|g" ../../Step_$(($Step+1))/Infos.dat
cp Best_Config_Ins.gro ../../Step_$(($Step+1))/Dynamic/$Project.gro
cp ../Dynamic/$Project.ndx ../../Step_$(($Step+1))/Dynamic
cp ../Dynamic/*.itp ../../Step_$(($Step+1))/Dynamic
#    cp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp_old
cp ../Dynamic/$Project.top ../../Step_$(($Step+1))/Dynamic
if [[ $updch == YES ]]; then
#    if [[ $MultChain == "YES" ]]; then
#       cp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp_old
#       cp chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp
#    else
#       cp ../../Step_$(($Step+1))/Dynamic/${Project}.top ../../Step_$(($Step+1))/Dynamic/${Project}.top_old
#       cp chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}.top
#    fi
    if [[ -f ../Dynamic/${Project}_Protein_chain_$chain.itp ]]; then
       cp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp_old
       cp chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}_Protein_chain_$chain.itp
    else
       cp ../../Step_$(($Step+1))/Dynamic/${Project}.top ../../Step_$(($Step+1))/Dynamic/${Project}.top_old
       cp chain.itp ../../Step_$(($Step+1))/Dynamic/${Project}.top
    fi
fi
cp -r ../Dynamic/amber94.ff ../../Step_$(($Step+1))/Dynamic/
cp ../Dynamic/residuetypes.dat ../../Step_$(($Step+1))/Dynamic
#cp $templatedir/standard-EM.mdp ../../Step_$(($Step+1))/Dynamic/
cp $templatedir/soglia ../../Step_$(($Step+1))/Dynamic/
#cp $templatedir/ASEC/DynIt_list_chromo.sh ../../Step_$(($Step+1))
cp ../Dynamic/$Project.ndx ../../Step_$(($Step+1))
cp ../$prm.prm ../../Step_$(($Step+1))
cp ../template_gro2tk ../../Step_$(($Step+1))
cp ../template_tk2gro ../../Step_$(($Step+1))

if [[ $solvent ==  "YES" ]]; then
   Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
   if [[ $moldy == "NPT" ]]; then
      cp $templatedir/ASEC/dynamic_sol_NPT.mdp ../../Step_$(($Step+1))/Dynamic/
   else
      cp $templatedir/ASEC/dynamic_sol_NVT.mdp ../../Step_$(($Step+1))/Dynamic/
   fi
else
   cp $templatedir/standard-EM.mdp ../../Step_$(($Step+1))/Dynamic/
fi
 
if [[ $Step -eq 0 ]]; then
   cp ../calculations/${Project}_$folder/${Project}_$folder.Final.xyz ../../Step_$(($Step+1))/${Project}_$folder.Final_last.xyz
   cp ../calculations/${Project}_$folder/${Project}_$folder.JobIph ../../Step_$(($Step+1))/${Project}_$folder.JobIph_old
   if [[ $uphess == YES ]]; then
      cp ../calculations/${Project}_$folder/${Project}_$folder.Hessian ../../Step_$(($Step+1))/${Project}_$folder.Hessian_old
   fi
else
   cp ../calculations/${Project}_$folder/${Project}_$folder.Final_last.xyz ../../Step_$(($Step+1))/${Project}_$folder.Final_last.xyz
   cp ../calculations/${Project}_$folder/${Project}_$folder.JobIph_new ../../Step_$(($Step+1))/${Project}_$folder.JobIph_old
   if [[ $uphess == YES ]]; then
      if [[ $fromhess == previous ]]; then
         cp ../calculations/${Project}_$folder/${Project}_$folder.Hessian_new_1 ../../Step_$(($Step+1))/${Project}_$folder.Hessian_old
      fi
      if [[ $fromhess == QMMM ]]; then
         cp ../${Project}_$folder.Hessian_QMMM ../../Step_$(($Step+1))/${Project}_$folder.Hessian_QMMM
      fi
   fi
fi

if [[ $solvent ==  "YES" ]]; then
   cp $templatedir/ASEC/MD_$moldy.sh ../../Step_$(($Step+1))
   echo ""    
   echo "Go to Step_$(($Step+1)) and continue with MD_$moldy.sh"
   echo ""
else
   cp $templatedir/ASEC/DynIt_list_chromo.sh ../../Step_$(($Step+1))
   echo ""    
   echo "Go to Step_$(($Step+1)) and continue with DinIt_list_chromo.sh"
   echo ""
fi

