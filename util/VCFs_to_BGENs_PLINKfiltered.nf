//Description : parallel conversion and filtering of VCF file to BGEN file

process Plinked {
  label "BGEN_generation"
  executor "local"
  cache "lenient"
    
  input :
  each file(VCF_file) from Channel.fromPath(params.VCF_files).flatten()
  each file(covar_file) from Channel.fromPath(params.covar_file)
  file pheno_file from Channel.fromPath(params.pheno_file)

  output:
  file "*.bgen" into bgen_file mode flatten
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
  """
name=${VCF_file.getSimpleName().replaceAll('.vcf$', '')}  

plink2 \
  --vcf $VCF_file \
  --pheno pheno_file \
  --covar covar_file \
  --maf ${params.maf} --geno ${params.geno} --hwe ${params.HWE} \
  --mind ${params.mind} \
  --max-alleles 2  \
  --export bgen-1.2 ${params.Plink2_Options} \
  --out \$name
  """
}
