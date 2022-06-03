#!/bin/bash


base_command="singularity run -B $InDir:$HOME/input  -B $OutDir:$HOME/ouput $SIF regenie \
    --step 1 \
    --loocv \
    --covarFile $HOME/input/$CovarName \
    --phenoFile $HOME/input/$PheName \
    --catCovarList $CovarCAT \
    --bsize $Bsize \
    --gz"
$base_command \
      --bed $HOME/input/$bfile \
      --out $HOME/ouput/test_bin_10 \
      --run-l0 $HOME/ouput/fit_bin_parallel.master,$SLURM_ARRAY_TASK_ID \
      --threads 4
