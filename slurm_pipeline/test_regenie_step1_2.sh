#!/bin/bash

base_command="singularity run -B $InDir:$HOME/input  -B $OutDir:$HOME/ouput $SIF regenie \
    --step 1 \
    --loocv \
    --covarFile $HOME/input/$CovarName \
    --phenoFile $HOME/input/$PheName \
    --bsize $Bsize \
    --gz"
    
$base_command \
    --bed $HOME/input/$bfile \
    --extract $HOME/ouput/"$qc_pass".snplist \
    --out $HOME/ouput/test_bin_10 \
    --run-l1 $HOME/ouput/fit_bin_parallel.master \
    --keep-l0 \
    --threads 5 \
    --force-step1
