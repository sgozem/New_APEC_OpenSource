#!/bin/bash
#
# Checking for the existence of Infos.dat
#
if [[ ! -f ../Infos.dat ]]; then
   echo " Fatal error! Infos.dat not found where I expected!"
   echo " Aborting..."
   echo " Minimization.sh 1 InfosNotFound" > ../arm.err
   exit 0
else
   echo " Minimization.sh 0 OK" > ../arm.err 
fi
#
# Checking for the Gromacs file required for minimizations
#
if [[ ! -f soglia ]]; then
   echo " Fatal error! soglia not found!"
   echo " Aborting..."
   echo ""
   echo " Minimization.sh 2 SogliaNotFound" > ../arm.err
   exit 0
fi
#
# Retrieving information from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`
gropath=`grep "GroPath" ../Infos.dat | awk '{ print $2 }'`
#
# Minimization is run until 2 consecutive steps has a RMSD difference less than
# the threshold specified in soglia
#
iter=1
flag=0
soglia=`cat soglia`
echo '0' > choices.txt
echo '0' >> choices.txt
while [ $flag -ne 1 ]; do
      $gropath/gmx mdrun -s $Project.tpr -o $Project.trr -x $Project.xtc -c final-$Project.gro 2> grolog
      control=`grep 'Can not' grolog`
      if [ ! -z $control ]; then
         echo " Fatal error! mdrun was not running properly!"
         echo ""
         echo " Minimization.sh 3 MdrunProblem" > ../arm.err
         exit 0
      fi
      $gropath/gmx g_rms -f $Project.trr -s $Project.tpr < choices.txt
      tail -n 50 rmsd.xvg > step50
      i=1
      while read line; do
            echo $line > linea
            rmsd[$i]=`awk '{ print $2 }' linea`
            if [ $i -gt 1 ]; then
               j=$(($i-1))
               diff=`echo ${rmsd[$i]}-${rmsd[$j]} | bc`
               if [[ "$diff" == *-* ]]; then
                  testo=`echo "-1*$diff < $soglia" | bc`
               else
                  testo=`echo "$diff < $soglia" | bc`
               fi
               if [[ $testo -eq 1 && $flag -ne 1 ]]; then
                  flag=1
                  steps=$i
               fi
            fi
            i=$(($i+1))
      done < step50
      rm linea step50
      if [ $flag -eq 1 ]; then
         echo " RMSD converged at step $steps of run $iter"
         echo ""
      else
         echo " RMSD not converged after 50 steps of run $iter"
         echo " mdrun will be executed for other 50 steps"
         echo ""
         sleep 1
         mkdir iter.$iter
         mv rmsd.xvg iter.$iter
         mv $Project.tpr iter.$iter
         mv $Project.trr iter.$iter
         mv ener.edr iter.$iter
         mv md.log iter.$iter
         mv grolog iter.$iter
         mv mdout.mdp iter.$iter
         mv $Project.gro iter.$iter
         cp final-$Project.gro iter.$iter
         mv final-$Project.gro $Project.gro
         $gropath/gmx grompp -f standard-EM.mdp -c $Project.gro -n $Project.ndx -p $Project.top -o $Project.tpr
      fi
      iter=$(($iter+1))
done
rm choices.txt

