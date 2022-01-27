#!/bin/bash

phe=$(head $InDir/$PheName -n 1 | tr ' ' '\t' | cut -f $SLURM_ARRAY_TASK_ID)


gunzip step2_*_$phe.regenie.gz
cat $OutDir/step2_1_"$phe".regenie > $OutDir/"$phe".regenie

for ((i=2; i<=$njobs; i++))
do
  sed 1d $OutDir/step2_"$i"_"$phe".regenie >> $OutDir/"$phe".regenie
done
