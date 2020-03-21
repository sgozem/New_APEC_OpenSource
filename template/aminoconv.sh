#!/bin/bash
#
lenout=$3
case $lenout in
     3) lett1=( $( echo $1 $2 ) )
        for ((i=0;$i<2;i=$(($i+1)))); do
            case ${lett1[$i]} in
                 G) lett3[$i]='GLY'
                 ;;
                 A) lett3[$i]='ALA'
                 ;;
                 V) lett3[$i]='VAL'
                 ;;
                 L) lett3[$i]='LEU'
                 ;;
                 I) lett3[$i]='ILE'
                 ;;
                 S) lett3[$i]='SER'
                 ;;
                 T) lett3[$i]='THR'
                 ;;
                 C) lett3[$i]='CYS'
                 ;;
                 M) lett3[$i]='MET'
                 ;;
                 P) lett3[$i]='PRO'
                 ;;
                 H) lett3[$i]='HIS'
                 ;;
                 R) lett3[$i]='ARG'
                 ;;
                 N) lett3[$i]='ASN'
                 ;;
                 Q) lett3[$i]='GLN'
                 ;;
                 E) lett3[$i]='GLU'
                 ;;
                 D) lett3[$i]='ASP'
                 ;;
                 F) lett3[$i]='PHE'
                 ;;
                 W) lett3[$i]='TRP'
                 ;;
                 Y) lett3[$i]='TYR'
                 ;;
                 K) lett3[$i]='LYS'
                 ;;
                 *) lett3[$i]='UNK'
                 ;;
             esac
        done
        echo ${lett3[0]} > wtres
        echo ${lett3[1]} > sidemut 
     ;;
     1) lett3=( $( echo $1 $2 ) )
        for ((i=0;$i<2;i=$(($i+1)))); do
            case ${lett3[$i]} in
                 GLY) lett1[$i]='G'
                 ;;
                 ALA) lett1[$i]='A'
                 ;;
                 VAL) lett1[$i]='V'
                 ;;
                 LEU) lett1[$i]='L'
                 ;;
                 ILE) lett1[$i]='I'
                 ;;
                 SER) lett1[$i]='S'
                 ;;
                 THR) lett1[$i]='T'
                 ;;
                 CYS) lett1[$i]='C'
                 ;;
                 MET) lett1[$i]='M'
                 ;;
                 PRO) lett1[$i]='P'
                 ;;
                 HIS) lett1[$i]='H'
                 ;;
                 ARG) lett1[$i]='R'
                 ;;
                 ASN) lett1[$i]='N'
                 ;;
                 GLN) lett1[$i]='Q'
                 ;;
                 GLU) lett1[$i]='E'
                 ;;
                 ASP) lett1[$i]='D'
                 ;;
                 PHE) lett1[$i]='F'
                 ;;
                 TRP) lett1[$i]='W'
                 ;;
                 TYR) lett1[$i]='Y'
                 ;;
                 LYS) lett1[$i]='K'
                 ;;
                 ASH) lett1[$i]='D'
                 ;;
                 GLH) lett1[$i]='E'
                 ;;
                 HID) lett1[$i]='H'
                 ;;
                 HIE) lett1[$i]='H'
                 ;;
                 LYN) lett1[$i]='K'
                 ;;
                 *) lett1[$i]='X'
                 ;;
             esac
        done
        echo ${lett1[0]} > wtres
        echo ${lett1[1]} > sidemut
     ;;
esac
