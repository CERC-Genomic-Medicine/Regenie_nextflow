//Description : parallel conversion and filtering of VCF file to BGEN file

process Plinked {
  label "BGEN_generation"
  cache "lenient"
  scratch true
    
  input :
  each file(VCF_file) from Channel.fromPath(params.VCF_files)


  output:
  file "*.bgen" into bgen_file mode flatten
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
  """
name=${VCF_file.getName().replaceAll('.vcf.gz$', '')}  

plink2 \
  --vcf $VCF_file \
  --maf ${params.maf} \ 
  --geno ${params.geno} \
  --hwe ${params.HWE} \
  --mind ${params.mind} \
  --max-alleles 2  \
  --export bgen-1.2 ${params.Plink2_Options} \
  --out \$name
  """
}
