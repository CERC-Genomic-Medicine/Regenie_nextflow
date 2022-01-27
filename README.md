# Regenie_Pipe
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

.regenie.gz file for each phenotype
