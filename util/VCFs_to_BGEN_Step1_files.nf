//Description : parallel conversion and filtering of VCF file to BGEN file

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

//Description : From vcf files (with indexes) to a single bgen file with custom filter

// Typicall VCF/BCF files are already split by chromosome. We make use of this, and process each of them in parallel.
process filter_by_chrom {
  cache "lenient"
  scratch true

  input:
  tuple file(vcf), file(vcf_index) from Channel.fromPath(params.VCF_files).map{ vcf -> [ vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".tbi") ] } // accepts VCF or BCF
  each file(Exclude) from excluded_unzip
  output:
  file "${vcf.getBaseName()}.common_independent_snps.vcf.*" into filtered_by_chrom mode flatten

  """
  # Apply hard filters and identify independent SNPs
  if [ ${vcf.getExtension()} = "bcf"  ];then
      plink2 \
      --bcf ${vcf} \
      --maf ${params.maf} \
      --geno ${params.geno} \
      --hwe ${params.HWE} \
      --min-alleles 2 \
      --max-alleles 2 \
      --snps-only \
      --set-all-var-ids '@_#_\$r_\$a' \
      --exclude ${Exclude} \
      --indep-pairwise 1000 100 0.9 \
      --export vcf-4.2 bgz ref-first \
      --out common_snps

  else
    plink2 \
      --vcf ${vcf} \
      --maf ${params.maf} \
      --geno ${params.geno} \
      --hwe ${params.HWE} \
      --min-alleles 2 \
      --max-alleles 2 \
      --snps-only \
      --set-all-var-ids '@_#_\$r_\$a' \
      --exclude ${Exclude} \
      --indep-pairwise 1000 100 0.9 \
      --export vcf-4.2 bgz ref-first \
      --out common_snps
  fi
  # Keep only independent SNPs in the final VCF
  plink2 \
    --vcf common_snps.vcf.gz \
    --extract common_snps.prune.in \
    --export vcf-4.2 bgz ref-first ${params.Plink2_Options}\
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
  file "all.common_independent_snps.bgen" into bgen
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
  file "all.common_independent_snps.sample" into sample
    publishDir "${params.OutDir}/BGEN", pattern: "*.sample", mode: "copy"
  file "all.common_independent_snps.bgen.bgi" into bgi_file mode flatten
    publishDir "${params.OutDir}/BGEN", pattern: "*.bgi", mode: "copy"
 script :   
  """
  find . -name "*.vcf.gz" | sort -V > files.txt
  bcftools concat -n -f files.txt -Oz -o all.common_independent_snps.vcf.gz
  gunzip -k all.common_independent_snps.vcf.gz
  bcftools index -t all.common_independent_snps.vcf.gz
  bash ${workflow.scriptFile.getParent()}/Scripts/dephasing.sh
  qctool_v2.2.0 -g all.common_independent_snps.vcf -filetype vcf -bgen-bits 8 -og all.common_independent_snps.bgen -os all.common_independent_snps.sample
  bgenix -index -g all.common_independent_snps.bgen
 """
}
