# Regenie_nextflow
Pipeline implementation of Regenie

## Before Starting ##
- Install this repository
```
git clone https://github.com/CERC-Genomic-Medicine/Regenie_pipe.git   ## the slurm_pipeline folder is depricated
```
- Install regenie image (if not already installed)  
```
module load singularity    
singularity pull docker://ghcr.io/rgcgithub/regenie/regenie:VERSION.gz    
tested with  singularity pull docker://ghcr.io/rgcgithub/regenie/regenie:v2.2.4.gz  
```

- Adjust the configuration file to your needs  

## Process ##  
```
nextflow run $PATH/Regenie.nf -c $PATH/nextflow.config

```

## Output ##

_[PHENOTYPE].regenie.gz file for each phenotype  
Chunk_*_.phe.txt for each phenotype parallelisation (phenotype bin)  
test_bin_[Phenotype_bin]_[Phenotype index in phenotype bin].loco.gz for each phenotype   
