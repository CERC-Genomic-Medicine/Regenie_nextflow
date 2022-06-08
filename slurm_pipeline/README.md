DEPRICATED

# Regenie_Pipe
Pipeline implementation of Regenie

## Beta version Caviat ##
Runs on slurm specifically Compute Canada due to dependency linked variables (see end of file for possible modification for further portability)  
Needs module Singularity plink V2 and recent boost.  
Includes a Plink filter (will be optionnal in futur versions)  
Future implementation will focus on scalability, portability, and the addition of more options  

## Before Starting ##
--install this repository (when public)  
```  
git clone https://github.com/CERC-Genomic-Medicine/Regenie_pipe.git   
```  

-Install regenie image (if not already installed)
```
module load singularity
```
download the singularity image of Regenie :
```
singularity pull docker://ghcr.io/rgcgithub/regenie/regenie:VERSION.gz
```
version tested with :  > singularity pull docker://ghcr.io/rgcgithub/regenie/regenie:v2.2.4.gz

-Adjust config file to your needs  
  -**Be carefull, any code added will run (will be fixed in futur versions)**  
 -Adjust time in Regine_pipeline_gen.sh (optional)

## Process ##
```
sbatch Regine_pipeline_gen.sh
```

## Output ##

.regenie file for each phenotype



## Portability caviat ##
To run on other slurm managed machine modify Regine_pipeline_gen.sh's jobnames variable to reflect the variable output of $(sbatch ...)
