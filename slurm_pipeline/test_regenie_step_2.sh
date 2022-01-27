#!/bin/bash

base_command="singularity run -B $InDir:$HOME/input  -B $OutDir:$HOME/ouput $SIF regenie \
    --step 2 \
    --lowmem \
    --covarFile $HOME/input/$CovarName \
    --phenoFile $HOME/input/$PheName \
    --bsize $Bsize \
    --gz"
    
$base_command \
    --bed $HOME/input/$bfile \
    --out $HOME/ouput/step2_$SLURM_ARRAY_TASK_ID \
    --pred $HOME/ouput/test_bin_10_pred.list \
    --extract $HOME/ouput/fit_bin_parallel_job$SLURM_ARRAY_TASK_ID.snplist \
    --threads 10
