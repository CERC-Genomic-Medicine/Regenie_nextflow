//Description : From vcf files (with indexes) to a single bgen file with custom filter

process Merge_Plinked {
  label "BGEN_generation_Merge_VCF"
  cache "lenient"
  scratch true
    
  input :
  file(VCF_file) from Channel.fromPath(params.VCF_files).collect()
  file(VCF_file_index) from Channel.fromPath(params.VCF_files_indexes).collect()
  val(out) from params.outname

  output:  
  file "*" into bgen_file mode flatten
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
    publishDir "${params.OutDir}/LOG", pattern: "*.log", mode: "copy"
    publishDir "${params.OutDir}/SAMPLE", pattern: "*.sample", mode: "copy"
    publishDir "${params.OutDir}/PRUNE/IN", pattern: "*.prune.in", mode: "copy" //optional
    publishDir "${params.OutDir}/PRUNE/OUT", pattern: "*.prune.out", mode: "copy" //optional
    
  """
bcftools concat ${VCF_file} -o VCF.vcf

plink2 \
  --vcf VCF.vcf \
  --maf ${params.maf} \
  --geno ${params.geno} \
  --mind ${params.mind} \
  --max-alleles 2  \
  --export bgen-1.2 ${params.Plink2_Options} \
  --out ${out}
  """
}
