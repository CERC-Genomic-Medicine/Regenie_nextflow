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
## Input file ##

P(B)GENs for variants to be tested (produced with util workflow from VCFs).  
 Sample file (if BGEN) or pvar & psam
P(B)GEN for common Variant (produced with util workflow from VCFs).   
 Sample file (if BGEN) or pvar & psam  
Covariant and Phenotype files (see format below)
```
FID IID Pheno1/Covar1 Pheno2/Covar2  
F1   S1      0.4            0.75
F2   S2      0.2             0.5
...
```
if there is no FID information use IIDx2 (and --double-id option in PGEN generation) psam must correspond to the FID ID column

## Process ##  
```
nextflow run $PATH/Regenie.nf -c $PATH/nextflow.config

```

# Regenie_nextflow/Util
useful pre-Regenie data handling  
software needed : BCFTOOLS, Plink version > 2, bgenix tool
  
  VCF(s)_to_B(P)GEN(s) usefull to produce pipeline input.
  
default config includes no filtering  
standard config is based on REGENIE paper (UK BioBank processing).

Usage
```
nextflow run $PATH/[script].nf -c $PATH/[config].config
```


