//Description : From vcf files (with indexes) to a single bgen file with custom filter

// Typicall VCF/BCF files are already split by chromosome. We make use of this, and process each of them in parallel.

process filter_by_chrom {
  cache "lenient"
  //scratch true

  input:
  file vcf from Channel.fromPath(params.VCF_files) // accepts VCF or BCF
  each file(exclude_bed) from Channel.fromPath("${workflow.projectDir}/Low_complexity_regions/${params.lcr_regions}.bed.gz")
  
  output:
  file "${vcf.getBaseName()}.common_independent_snps.*" into filtered_by_chrom mode flatten
      
  """
  # Apply hard filters and identify independent SNPs
  if [ ${vcf.getExtension()} = "bcf" ]; then
     plink_import_option="--bcf ${vcf}"
  else
     plink_import_option="--vcf ${vcf}"
  fi

  ${params.plink2_exec} \${plink_import_option} \
    --maf ${params.maf} \
    --geno ${params.geno} \
    --hwe ${params.HWE} \
    --min-alleles 2 \
    --max-alleles 2 \
    --exclude bed0 ${exclude_bed} \
    --snps-only \
    --set-all-var-ids '@:#:\$r:\$a' \
    --indep-pairwise 1000 100 0.9 \
    --make-pgen \
    --out common_snps

  # Keep only independent SNPs
  ${params.plink2_exec} \
    --pfile common_snps \
    --extract common_snps.prune.in \
    --make-pgen erase-phase \
    --out ${vcf.getBaseName()}.common_independent_snps
  """
}


process merge_chroms {
  cache "lenient"
  //scratch true

  input:
  file(pfiles) from filtered_by_chrom.collect()

  output:
  file "all.common_independent_snps.*" into pgen

  publishDir "${params.OutDir}/", pattern: "all.common_independent_snps.*", mode: "copy"
    
  script:
  if (params.format == "PGEN")
     """
     find . -name "*.pgen" -printf "%f\n" | sort -V | sed s"/.pgen//" > files.txt
     ${params.plink2_exec} --pmerge-list files.txt --make-pgen --out all.common_independent_snps
     """

  else if (params.format == "BGEN")
     """
     find . -name "*.pgen" -printf "%f\n" | sort -V | sed s"/.pgen//" > files.txt
     ${params.plink2_exec} --pmerge-list files.txt --export vcf-4.2 bgz ref-first --out temporary_merged_vcf
     ${params.qctool_exec} -g temporary_merged_vcf.vcf.gz -filetype vcf -bgen-bits 8 -og all.common_independent_snps.bgen -os all.common_independent_snps.sample
     ${params.bgenix_exec} -index -g all.common_independent_snps.bgen
     """
  
  else
     error "Invalid format: ${params.format}"
}
