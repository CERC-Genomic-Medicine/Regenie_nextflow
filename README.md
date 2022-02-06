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

# Regenie_nextflow/Util
useful pre-Regenie data handling
software needed : BCFTOOLS and Plink version > 2  

 - VCFs_to_BGEN_PLINKfiltered.nf             //VCFs transition to BGEN (many to one) with filters [config default no filter]  
 - VCFs_to_BGENs_PLINKfiltered.nf            //VCFs transition to BGENs (many to many) with filters [config default no filter]

modify util.config for filtering

Usage
```
nextflow run $PATH/[script].nf -c $PATH/util.config
```


