#!/bin/bash

phe=$(head $InDir/$PheName -n 1 | tr ' ' '\t' | cut -f $SLURM_ARRAY_TASK_ID)

gunzip step2_*_$phe.regenie.gz
cat step2_1_$phe.regenie > $phe.regenie

for ((i=2; i<=$njobs; i++))
do
  sed 1d step2_$i_$phe.regenie >> $phe.regenie
done
