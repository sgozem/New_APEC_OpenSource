#!/bin/bash
#
# @ shell       = /bin/bash
# @ job_name    = NOMEPROGETTO
# @ output      = NOMEPROGETTO.out
# @ error       = NOMEPROGETTO.err
# @ wall_clock_limit = 144:00:00
# @ job_type    = serial
# @ resources = ConsumableMemory(800Mb)
# NB: The following line for ISCRA account holders only
# @ account_no = IscrA_ASREXMD
# @ queue
export inpdir=NOMEDIRETTORI
export outdir=NOMEDIRETTORI/output
cp $inpdir/* $WorkDir
cd $WorkDir
GROPATH/mdrun -s $Project.tpr -o $Project.trr -x $Project.xtc -c final-$Project.gro
cp * $outdir/

