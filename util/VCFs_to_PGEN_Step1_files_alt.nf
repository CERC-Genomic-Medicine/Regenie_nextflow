//Description : From vcf files (with indexes) to a single bgen file with custom filter
process Unizping {

input:
file(exclude) from Channel.fromPath(params.bed)

output:
  file("*.bed") into excluded_unzip

  """
if [${params.bed[-6..-1]}=="tar.gz"]; then
tar -xzvf ${exclude}
else
cp ${exclude} file.bed
fi
  """
}

// Typicall VCF/BCF files are already split by chromosome. We make use of this, and process each of them in parallel.
process filter_by_chrom {
  cache "lenient"
  scratch true

  input:
  tuple file(vcf), file(vcf_index) from Channel.fromPath(params.VCF_files).map{ vcf -> [ vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".tbi") ] } // accepts VCF or BCF
  each file(Exclude) from excluded_unzip
  
  output:
  file "${vcf.getBaseName()}.common_independent_snps.vcf.*" into filtered_by_chrom mode flatten
    publishDir "${params.OutDir}/VCF", pattern: "*.vcf.gz*", mode: "copy"
  """
  # Apply hard filters and identify independent SNPs
  if [ ${vcf.getExtension()} = "bcf" ];then
      plink2 \
      --bcf ${vcf} \
      --maf ${params.maf} \
      --geno ${params.geno} \
      --hwe ${params.HWE} \
      --min-alleles 2 \
      --max-alleles 2 \
      --snps-only \
      --set-all-var-ids '@_#_\$r_\$a' \
      --indep-pairwise 1000 100 0.9 \
      --export vcf-4.2 bgz ref-first \
      --exclude ${Exclude} \
      --out common_snps

  else 
    plink2 \
      --vcf ${vcf} \
      --maf ${params.maf} \
      --geno ${params.geno} \
      --hwe ${params.HWE} \
      --min-alleles 2 \
      --max-alleles 2 \
      --exclude ${Exclude} \
      --snps-only \
      --set-all-var-ids '@_#_\$r_\$a' \
      --indep-pairwise 1000 100 0.9 \
      --export vcf-4.2 bgz ref-first \
      --out common_snps
  fi
  # Keep only independent SNPs in the final VCF
  plink2 \
    --vcf common_snps.vcf.gz \
    --extract common_snps.prune.in \
    --export vcf-4.2 bgz ref-first \
    --out ${vcf.getBaseName()}.common_independent_snps
  
  # Create index. It will be used when concatinating.
  bcftools index -t ${vcf.getBaseName()}.common_independent_snps.vcf.gz
  """
}


process merge {
  cache "lenient"
  scratch true

  input:
  file(vcfs_and_indices) from filtered_by_chrom.collect()

  output:
  tuple file("all.common_independent_snps.vcf.gz"), file("all.common_independent_snps.vcf.gz.tbi") into vcf_final mode flatten
    publishDir "${params.OutDir}/VCF", pattern: "*.vcf.gz*", mode: "copy"
  file "*.pgen" into pgen_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.pgen", mode: "copy"
  file "*.pvar" into pvar_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.pvar", mode: "copy"
  file "*.psam" into psam_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.psam", mode: "copy"
    
 script :   
  """
  find . -name "*.vcf.gz" | sort -V > files.txt
  bcftools concat -n -f files.txt -Oz -o all.common_independent_snps.vcf.gz
  tabix -t all.common_independent_snps.vcf.gz
  plink2 --vcf all.common_independent_snps.vcf.gz --make-pgen 'erase-phase' ${params.Plink2_Options} --out all.common_independent_snps
 """
}
