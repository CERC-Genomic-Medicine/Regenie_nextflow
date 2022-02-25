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
software needed : QCTOOOL, Plink version > 2

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
or if you want to use the pgen version (requieres pgen input files)
```
nextflow run $PATH/Regenie_vpgen.nf -c $PATH/nextflow.config

```

# Regenie_nextflow/Util
useful pre-Regenie data handling  
software needed : QCTOOOL, Plink version > 2, BCFTOOLS
  
  VCF(s)_to_B(P)GEN(s) usefull to produce pipeline input.
  
default config includes no filtering  
standard config is based on REGENIE paper (UK BioBank processing).

Usage
```
nextflow run $PATH/[script].nf -c $PATH/[config].config
```

# Necessary software and their installation

BGEN most recent version available @ https://enkre.net/cgi-bin/code/bgen/dir?ci=tip (module load nixpkgs/16.09 gcc/5.4.0 bgen/1.1.4)
Plink2 avaullable @ https://github.com/chrchang/plink-ng
QCTOOL 
