#! /bin/bash

source parameters.config

# first job - no dependencies
jid1=$(sbatch --account=$nameaccount --job-name=regenie_step_0 --mem=10G --cpus-per-task=1 --time=10:15:00 $wkdir/test_regenie_step1_0.sh)
job1=$(echo $jid1 | cut -f 4 -d " " -) ### Compute canada slurm does "Submitted batch job jobID" we only want jobID

# multiple jobs can depend on a single job
jid2=$(sbatch --account=$nameaccount --dependency=afterok:$job1 --cpus-per-task=4 --array=1-$njobs%10 --time=10:15:00 --mem-per-cpu=2g $wkdir/test_regenie_step1_1.sh)
job2=$(echo $jid2 | cut -f 4 -d " " -) ### Compute canada slurm does "Submitted batch job jobID" we only want jobID

# multiple jobs can depend on a single job
jid3=$(sbatch --account=$nameaccount --dependency=afterok:$job2 --cpus-per-task=5 --time=10:15:00 --mem-per-cpu=8g $wkdir/test_regenie_step1_2.sh)
job3=$(echo $jid3 | cut -f 4 -d " " -) ### Compute canada slurm does "Submitted batch job jobID" we only want jobID

# multiple jobs can depend on a single job
jid4=$(sbatch --account=$nameaccount --dependency=afterok:$job3 --cpus-per-task=10 --array=1-$njobs%10 --time=10:15:00 --mem-per-cpu=8g $wkdir/test_regenie_step_2.sh)
job4=$(echo $jid4 | cut -f 4 -d " " -) ### Compute canada slurm does "Submitted batch job jobID" we only want jobID

#Concat
PHENO=$(head $InDir/$PheName -n 1 |  tr '\t' ' ' | wc -w)
jid5=$(sbatch --account=$nameaccount --dependency=afterok:$job4 --cpus-per-task=10 --array=3-$PHENO%10 --time=10:15:00 --mem-per-cpu=8g $wkdir/test_regenie_concat.sh)
