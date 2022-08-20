#!/bin/bash
#
# This script retrieves CASPT2 energies from CASPT2 single point output file
# First of all, the project name is retrieved from Infos.dat
#
Project=`grep "Project" ../Infos.dat | awk '{ print $2 }'`

# Going to CASPT2 and grepping the data from the output file
#
caspt2=${Project}_CASPT2
cd $caspt2
grep -A3 "  Total CASPT2 energies:" $caspt2.out > ../CASPT2_Energies
grep -A10 "Dipole transition strengths" $caspt2.out | grep -B2 '2    3' | awk '{ print $3 }' > ../f_Oscillator
cd ..

# Getting lines 2, 3 and 4 of CASPT2_Energies
#
qw=CASPT2_Energies

s0=`sed -n 2p $qw | awk '{ print $7 }'`
s1=`sed -n 3p $qw | awk '{ print $7 }'`
s2=`sed -n 4p $qw | awk '{ print $7 }'`

# Output of the energies (a.u.) in the CASPT2_Energies file
#
echo "" >> $qw
echo " CASPT2 energies:" >> $qw
echo " S0 = $s0 a.u." >> $qw
echo " S1 = $s1 a.u." >> $qw
echo " S2 = $s2 a.u." >> $qw

# Calculating the energy gaps and converting them in wavelengths and kcal mol-1
# The S0-S1 and S0-S2 are calculated and printed
#
delta0_1=`echo "$s1 - $s0" | bc`
delta0_2=`echo "$s2 - $s0" | bc`
lambda0_1=`echo "45.5628759223273 / $delta0_1" | bc`
lambda0_2=`echo "45.5628759223273 / $delta0_2" | bc`
kcal0_1=`echo "627.50947 * $delta0_1" | bc`
kcal0_2=`echo "627.50947 * $delta0_2" | bc`

# Output of the calculated gaps in CASPT2_Energies
#
echo "" >> $qw
echo " Maximum absorption wavelength" >> $qw
echo " Root1_2 = $lambda0_1 nm" >> $qw
echo " Root1_3 = $lambda0_2 nm" >> $qw
echo "" >> $qw
echo " Energy gap" >> $qw
echo " Root1_2 = $kcal0_1 kcal/mol" >> $qw
echo " Root1_3 = $kcal0_2 kcal/mol" >> $qw
echo "" >> $qw

# Gettin oscillator strengths from f_Oscillator
#
oscil=f_Oscillator

onetwo=`sed -n 1p $oscil`
onethree=`sed -n 2p $oscil`
twothree=`sed -n 3p $oscil`

# Output of oscillator strengths in f_Oscillator
#
rm $oscil
echo "" >> $oscil
echo " Oscillator strengths at the CASSCF level:" >> $oscil
echo " Root1 -> Root2 = $onetwo" >> $oscil
echo " Root1 -> Root3 = $onethree" >> $oscil
echo " Root2 -> Root3 = $twothree" >> $oscil
echo "" >> $oscil

# Messages to the user
#
echo ""
echo " All the CASPT2 results were saved to $qw"
echo " The oscillator strengths were saved to $oscil"
echo ""

