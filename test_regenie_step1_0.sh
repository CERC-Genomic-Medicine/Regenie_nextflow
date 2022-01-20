#!/bin/bash


plink2 \
  --bfile $InDir/$bfile \
  --maf $maf --mac $mac --geno $geno --hwe $HWE \
  --mind $mind \
  --write-snplist --write-samples --no-id-header \
  --out $OutDir/$qc_pass


base_command="singularity run -B $InDir:$HOME/input  -B $OutDir:$HOME/ouput $SIF regenie \
    --step 1 \
    --loocv \
    --covarFile $HOME/input/$CovarName \
    --phenoFile $HOME/input/$PheName \
    --bsize $Bsize \
    --gz"
    
$base_command \
    --bed $HOME/input/$bfile \
    --out $HOME/ouput/test_bin_10 \
    --split-l0 $HOME/ouput/fit_bin_parallel,$njobs \
    --extract $HOME/ouput/"$qc_pass".snplist \
    --threads 1 \
    --force-step1
