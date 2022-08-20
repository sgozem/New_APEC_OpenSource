#!/bin/bash
#
# @ shell       = /bin/bash
# @ job_name    = NOMEPROGETTO
# @ output      = NOMEPROGETTO.$(jobid).out
# @ error       = NOMEPROGETTO.$(jobid).err
# @ wall_clock_limit = hh:00:00
# @ job_type    = serial
# @ resources = ConsumableMemory(MEMTOTMb)
# NB: The following line for ISCRA account holders only
# @ account_no = IscrA_ASREXMD
# @ queue

export Project=$LOADL_JOB_NAME
export MOLCAS=/nfs/03/bgs0361/bin/7.8.dev
export MOLCASMEM=MEMORIAMB
export TINKER=/nfs/03/bgs0361/bin/7.8.dev/tinker/bin_qmmm
export WorkDir=$CINECA_SCRATCH/$Project.$RANDOM
export InpDir=NOMEDIRETTORI

/nfs/03/bgs0361/bin/dowser/bin/molcas $InpDir/$Project.input > $InpDir/$Project.out 2> $InpDir/$Project.err

