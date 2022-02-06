//Description : From vcf files (with indexes) to a single bgen file with custom filter

process Merge_Plinked {
  label "BGEN_generation_Merge_VCF"
  executor "local"
  cache "lenient"
    
  input :
  file(VCF_file) from Channel.fromPath(params.VCF_files).flatten().collect()
  file(VCF_file_index) from Channel.fromPath(params.VCF_files_indexes).flatten().collect()
  file pheno_file from Channel.fromPath(params.pheno_file)
  val(out) from params.outname

  output:
  file "*.bgen" into bgen_file mode flatten
  publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
  """
bcftools concat ${VCF_file} -o VCF.vcf

plink2 \
  --vcf VCF.vcf \
  --pheno ${pheno_file} \
  --maf ${params.maf} --geno ${params.geno} --hwe ${params.HWE} \
  --mind ${params.mind} \
  --max-alleles 2  \
  --export bgen-1.2 ${params.Plink2_Options} \
  --out ${out}
  """
}
