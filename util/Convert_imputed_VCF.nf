//Description : parallel conversion and filtering of VCF file to PGEN file

process convert {
  cache "lenient"
  //scratch true
  
  input :
  tuple val(name), path(vcf)

  output:
  path("${name}.${params.vcf_field}.*")

  publishDir "${params.OutDir}/", pattern: "${name}.${params.vcf_field}.*", mode: "copy"
 
  script:
  if (params.format == "PGEN")
     """
     if [ ${vcf.getExtension()} = "bcf" ]; then
        plink_import_option="--bcf ${vcf}"
     else
        plink_import_option="--vcf ${vcf}"
     fi

     if [ ${params.vcf_field} != "GT" ]; then
        plink_import_option="\${plink_import_option} dosage=${params.vcf_field}"
     fi
     
     ${params.plink2_exec} \${plink_import_option} ${params.name} --max-alleles 2 --new-id-max-allele-len 8000 --set-all-var-ids '@:#:\$r:\$a' --make-pgen erase-phase --out ${name}.${params.vcf_field}
     """

  else if (params.format == "BGEN")
     """
     qctool_export_option=""
     if [ ${params.vcf_field} != "GT" ]; then
        qctool_export_option="-vcf-genotype-field ${params.vcf_field}"
     fi

     ${params.bcftools_exec} view -m2 -M2 ${vcf} | ${params.bcftools_exec} annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' -x ^FORMAT/GT,FORMAT/GP | sed -e s"/\\(\\s[[:digit:]]\\|\\s\\.\\)|\\([[:digit:]]:\\|\\.:\\)/\\1\\/\\2/g" | ${params.qctool_exec} -g - -filetype vcf \${qctool_export_option} -bgen-bits 8 -og ${name}.${params.vcf_field}.bgen -os ${name}.${params.vcf_field}.sample
     ${params.bgenix_exec} -index -g ${name}.${params.vcf_field}.bgen
     """

  else
     error "Invalid format: ${params.format}"
}


workflow {
input_vcf=Channel.fromPath(params.VCF_files).map { vcf -> [vcf.getName().replaceAll('.vcf.gz$|.bcf$', ''), vcf ] } // accepts VCF or BCF

convert(input_vcf)

}
