#!/bin/bash

phe=$(head $InDir/$PheName -n 1 | tr ' ' '\t' | cut -f $SLURM_ARRAY_TASK_ID)
echo $phe
echo $SLURM_ARRAY_TASK_ID
gunzip step2_*_$phe.regenie.gz
cat step2_*_$phe.regenie > $phe.regenie
