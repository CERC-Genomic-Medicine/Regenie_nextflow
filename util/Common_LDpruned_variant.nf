#!/usr/bin/env nextflow

/*
* AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2024
*/

// This pipeline performs filtering and  LD pruning of genetic variants to be used for PCA or for computing genomic predictions in Regenie.
// The input files are expected to be already split by chromosome for parallel processing.
// If input genotypes are in VCF/BCF, then males (if any) in chromosome X non-PAR regions must be coded as haploids (i.e. 0 and 1, instead of 0/0 and 1/1).


process filter_by_chrom_VCF_BCF {
   cache "lenient"
   //scratch true

   cpus 4
   memory "16 GB"
   time "1h"

   input:
   tuple val(filename), path(accompanying_files)
  
   output:
   path "${filename}.common_independent_snps.p*",  emit: out
   path "*.log", emit: logs

   publishDir "${params.output_dir}/logs", pattern: "*.log", mode: "copy"

   """
   # Set import option
   if [[ ${filename} =~ "bcf\$" ]]; then
      plink_import_option="--bcf ${filename} ${params.name}"
   else
      plink_import_option="--vcf ${filename} ${params.name}"
   fi
  
   # Exclude long-range LD and low complexity regions
   plink_exclude_option="--exclude bed0 ${workflow.projectDir}/Low_complexity_regions/${params.ld_regions}.bed.gz ${workflow.projectDir}/Low_complexity_regions/${params.lcr_regions}.bed.gz"
  
   if [[ -z "${params.hapmap_sites}" ]]; then
      plink_extract_option=""
   else
      plink_extract_option="--extract bed0 ${workflow.projectDir}/Low_complexity_regions/${params.hapmap_sites}.bed.gz"
   fi

   # Check if the file contains chromosome X non-PAR region. If it does, then extract ploidy information to pass to PLINK2.
   has_nonpar=`bcftools view -HG ${filename} ${params.nonpar_region} | head -n1 | wc -l`
   if [[ "\$has_nonpar" -eq 1 ]]; then
      bcftools view ${filename} ${params.nonpar_region} | bcftools +check-ploidy | tail -n+2 > ploidy.txt
      
      # Sanity check: individual IDs (their order and number)  in ploidy file must match individuals IDs in VCF header. If not, some individuals have mixed ploidy and were reported multiple times in the ploidy file.
      cut -f1 ploidy.txt > samples_ploidy.txt
      bcftools query -l ${filename} > samples_vcf.txt
      if ! cmp -s samples_ploidy.txt samples_vcf.txt; then
         exit 1 # there are individuals with different ploidy at different variants
      fi

      echo -e "#FID\tIID\tSEX" > ploidy.psam
      awk '{OFS = "\t" ; print \$1, \$1, \$5}' ploidy.txt >> ploidy.psam 
      
      plink_chrX_options="--psam ploidy.psam --split-par hg38"
   else
      plink_chrX_options=""
   fi

  # Apply filters and identify independent SNPs
  ${params.plink2_exec} \${plink_import_option} \
    \${plink_exclude_option} \
    \${plink_extract_option} \
    \${plink_chrX_options} \
    --geno ${params.geno} \
    --maf ${params.maf} \
    --min-alleles 2 \
    --max-alleles 2 \
    --snps-only \
    --hwe ${params.HWE} \
    --set-all-var-ids '@:#:\$r:\$a' \
    --indep-pairwise 1000 100 ${params.Rsq} \
    --make-pgen \
    --out ${filename}.common_snps

  # Keep only independent SNPs
  ${params.plink2_exec} \
    --pfile ${filename}.common_snps \
    --extract ${filename}.common_snps.prune.in \
    --make-pgen erase-phase \
    --out ${filename}.common_independent_snps
  """
}


process filter_by_chrom_pgen_bgen {
   cache "lenient"
   //scratch true

   cpus 4
   memory "16 GB"
   time "1h"

   input:
   tuple val(prefix), path(accompanying_files)
  
   output:
   path "${prefix}.common_independent_snps.p*", emit: out
   path "*.log", emit: logs

   """
   # Set import option
   if [ -f "${prefix}.pgen" ]; then
      plink_import_option="--pgen ${prefix}"
   else
      plink_import_option="--bgen ${prefix}"
   fi

   # Exclude long-range LD and low complexity regions
   plink_exclude_option="--exclude bed0 ${workflow.projectDir}/Low_complexity_regions/${params.ld_regions}.bed.gz ${workflow.projectDir}/Low_complexity_regions/${params.lcr_regions}.bed.gz"

   if [[ -z "${params.hapmap_sites}" ]]; then
      plink_extract_option=""
   else
      plink_extract_option="--extract bed0 ${workflow.projectDir}/Low_complexity_regions/${params.hapmap_sites}.bed.gz"
   fi

   # Apply filters and identify independent SNPs.
   # Note that here we do not do anything special about chromosome X, since provided .psam and .sample files allow PLINK2 to handle chromosome X ploidy automatically. Remember that it is up to you to confirm that your input .psam and .sample files have correct sex info.
   ${params.plink2_exec} \${plink_import_option} \
     \${plink_exclude_option} \
     \${plink_extract_option} \
     --geno ${params.geno} \
     --maf ${params.maf} \
     --min-alleles 2 \
     --max-alleles 2 \
     --snps-only \
     --hwe ${params.HWE} \
     --set-all-var-ids '@:#:\$r:\$a' \
     --indep-pairwise 1000 100 ${params.Rsq} \
     --make-pgen \
     --out ${prefix}.common_snps

   # Keep only independent SNPs
   ${params.plink2_exec} \
     --pfile ${prefix}.common_snps \
     --extract ${prefix}.common_snps.prune.in \
     --make-pgen erase-phase \
     --out ${prefix}.common_independent_snps
  """
}


process merge_chroms {
   cache "lenient"
   //scratch true

   cpus 1
   memory "4 GB"
   time "1h"

   input:
   path(pfiles)

   output:
   path "all.common_independent_snps.*"
   
   publishDir "${params.output_dir}/", pattern: "all.common_independent_snps.*", mode: "copy"
    
   script:
   if (params.output_format == "PGEN")
      """
      find . -name "*.pgen" -printf "%f\n" | sort -V | sed s"/.pgen//" > files.txt
      ${params.plink2_exec} --pmerge-list files.txt --make-pgen ${params.name} --out all.common_independent_snps
      """
   else if (params.output_format == "BGEN")
      """
      find . -name "*.pgen" -printf "%f\n" | sort -V | sed s"/.pgen//" > files.txt
      ${params.plink2_exec} --pmerge-list files.txt --export vcf-4.2 bgz ref-first --out temporary_merged_vcf
      ${params.qctool_exec} -g temporary_merged_vcf.vcf.gz -filetype vcf -bgen-bits 8 -og all.common_independent_snps.bgen -os all.common_independent_snps.sample
      ${params.bgenix_exec} -index -g all.common_independent_snps.bgen
      """
   else
      error "Invalid output format: ${params.output_format}"
}


workflow {
   // This channel will hold the following tuples:
   // If VCF/BCF: [filename, [/full/path/to/filename, /full/path/to/filename.tbi]] or [filename, [/full/path/to/filename, /full/path/to/filename.csi]]
   // If PGEN: [prefix, [/full/path/to/prefix.pgen, /full/path/to/prefix.psam, /full/path/to/prefix.pvar]]
   // If BGEN: [prefix, [/full/path/to/prefix.bgen, /full/path/to/prefix.sample]]
   genotype_files = Channel.fromFilePairs("${params.genotype_files}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi|.psam|.pvar|.pgen|.bgen|.sample)\$", "")}
   
   // Separate VCF/BCF from PLINK/BGEN
   genotype_files.branch {
      vcfs: it[0] =~ /.vcf.gz|bcf/
      non_vcfs: true
   }.set{ genotype_files_by_type }
   //genotype_files_by_type.vcf.view()


   // Instead of if statement lets use the concat. Assumption is that only one of the channels is non-empty (i.e. mixture of different file types is highly unlikely and should not happen).
   filtered_by_chrom = filter_by_chrom_VCF_BCF(genotype_files_by_type.vcfs).out.concat(filter_by_chrom_pgen_bgen(genotype_files_by_type.non_vcfs).out)
   //filtered_by_chrom.view()

   merge_chroms(filtered_by_chrom.collect())
}
