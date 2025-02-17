#!/usr/bin/env nextflow

/*
* AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2025
*/


// This pipeline converts VCF/BCF files imputed using minimac (or TOPMed/Michigan imputation servers) to the PGEN (PLINK2) or BGEN formats.


process convert {
  cache "lenient"
  //scratch true

   cpus 4
   memory "64 GB"
   time "48h"
  
  input :
  tuple val(filename), path(accompanying_files)

  output:
  path("${filename}.${params.vcf_field}.*")

  publishDir "${params.output_dir}", pattern: "${filename}.${params.vcf_field}.*"
 
  script:
  if (params.output_format == "PGEN")
     """
     if [[ ${filename} =~ "bcf\$" ]]; then
        plink_import_option="--bcf ${filename}"
     else
        plink_import_option="--vcf ${filename}"
     fi

     if [ ${params.vcf_field} != "GT" ]; then
        plink_import_option="\${plink_import_option} dosage=${params.vcf_field}"
     fi

     # Check if the file contains chromosome X non-PAR region. If it does, then extract ploidy information to pass to PLINK2.
     has_nonpar=`${params.bcftools_exec} view -HG ${filename} ${params.nonpar_region} | head -n1 | wc -l`
     if [[ "\$has_nonpar" -eq 1 ]]; then
        ${params.bcftools_exec} annotate -xINFO,FORMAT -Ou -r ${params.nonpar_region} ${filename} | bcftools +check-ploidy | tail -n+2 > ploidy.txt
	
        # Sanity check: individual IDs (their order and number)  in ploidy file must match individuals IDs in VCF header. If not, some individuals have mixed ploidy and were reported multiple times in the ploidy file.
        cut -f1 ploidy.txt > samples_ploidy.txt
        ${params.bcftools_exec} query -l ${filename} > samples_vcf.txt
        if ! cmp -s samples_ploidy.txt samples_vcf.txt; then
           exit 1 # there are individuals with different ploidy at different variants
        fi

        echo -e "#FID\tIID\tSEX" > ploidy.psam
        awk '{OFS = "\t" ; print \$1, \$1, \$5}' ploidy.txt >> ploidy.psam

	plink_chrX_options="--psam ploidy.psam --split-par hg38"
     else
        plink_chrX_options=""
     fi

     ${params.plink2_exec} \${plink_import_option} \${plink_chrX_options} ${params.name} --max-alleles 2 --new-id-max-allele-len 8000 --set-all-var-ids '@:#:\$r:\$a' --make-pgen erase-phase --out ${filename}.${params.vcf_field}
     """

  else if (params.output_format == "BGEN")
     """
     qctool_export_option=""
     if [ ${params.vcf_field} != "GT" ]; then
        qctool_export_option="-vcf-genotype-field ${params.vcf_field}"
     fi

     ${params.bcftools_exec} view -m2 -M2 -Ou ${filename} | ${params.bcftools_exec} annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' -x ^FORMAT/GT,FORMAT/GP | sed -e s"/\\(\\s[[:digit:]]\\|\\s\\.\\)|\\([[:digit:]]:\\|\\.:\\)/\\1\\/\\2/g" | ${params.qctool_exec} -g - -filetype vcf \${qctool_export_option} -bgen-bits 8 -og ${name}.${params.vcf_field}.bgen -os ${name}.${params.vcf_field}.sample
     ${params.bgenix_exec} -index -g ${filename}.${params.vcf_field}.bgen
     """

  else
     error "Invalid format: ${params.output_format}"
}


workflow {
   // This channel will hold the following tuples:
   // If VCF/BCF: [filename, [/full/path/to/filename, /full/path/to/filename.tbi]] or [filename, [/full/path/to/filename, /full/path/to/filename.csi]]
   input_files = Channel.fromFilePairs("${params.imputed_genotype_files}", size: -1) { file -> file.getName().replaceAll("(.tbi|.csi)\$", "")}

   //input_files.view()

   convert(input_files)
}
