//Description : parallel conversion and filtering of VCF file to BGEN file

process Plinked {
  label "BGEN_generation"
  cache "lenient"
  scratch true
    
  input :
  each file(VCF_file) from Channel.fromPath(params.VCF_files)


  output:
  output:  
  file "*" into bgen_file mode flatten
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
    publishDir "${params.OutDir}/LOG", pattern: "*.log", mode: "copy"
    publishDir "${params.OutDir}/SAMPLE", pattern: "*.sample", mode: "copy"
    publishDir "${params.OutDir}/PRUNE/IN", pattern: "*.prune.in", mode: "copy" //optional
    publishDir "${params.OutDir}/PRUNE/OUT", pattern: "*.prune.out", mode: "copy" //optional
  """
name=${VCF_file.getName().replaceAll('.vcf.gz$', '')}  

plink2 \
  --vcf $VCF_file \
  --maf ${params.maf} \ 
  --geno ${params.geno} \
  --mind ${params.mind} \
  --max-alleles 2  \
  --export bgen-1.2 ${params.Plink2_Options} \
  --out \$name
  """
}
