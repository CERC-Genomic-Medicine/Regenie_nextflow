//Description : parallel conversion and filtering of VCF file to PGEN file


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

process Plinked {
  label "PGEN_generation"
  cache "lenient"
  scratch true
    
  input :
  tuple file(vcf), file(vcf_index) from Channel.fromPath(params.VCF_files).map{ vcf -> [ vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".tbi") ] } // accepts VCF or BCF
  each file(Exclude) from excluded_unzip

  output:
  file "*.pgen" into pgen_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.pgen", mode: "copy"
  file "*.pvar" into pvar_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.pvar", mode: "copy"
  file "*.psam" into psam_file mode flatten
      publishDir "${params.OutDir}/pfile", pattern: "*.psam", mode: "copy"
    
  """
  if [ ${vcf.getExtension()} = "bcf"  ] ; then
    name=${vcf.getName().replaceAll('.bcf$', '')}
   if [ ${params.dosageFD} = "none"  ];then
    plink2 --bcf ${vcf}  --max-alleles 2 --exclude ${Exclude}  --new-id-max-allele-len 1000 --set-all-var-ids '@_#_\$r_\$a' --make-pgen 'erase-phase' ${params.Plink2_Options} --out \$name
   else
    plink2 --bcf ${vcf}  'dosage=${params.dosageFD}' --exclude ${Exclude} --new-id-max-allele-len 1000  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' ${params.Plink2_Options} --make-pgen 'erase-phase' --out \$name
   fi
  else
   if [ ${vcf.getExtension()} = "gz"  ];then
    name=${vcf.getName().replaceAll('.vcf.gz$', '')}
   else
    name=${vcf.getName().replaceAll('.vcf$', '')}  
   fi  
   if [ ${params.dosageFD} = "none"  ];then
    plink2 --vcf ${vcf}  --max-alleles 2 --exclude ${Exclude} --make-pgen 'erase-phase'  --set-all-var-ids '@_#_\$r_\$a' ${params.Plink2_Options} --new-id-max-allele-len 1000 --out \$name
   else
    plink2 --vcf ${vcf}  'dosage=${params.dosageFD}' --exclude ${Exclude}  --new-id-max-allele-len 1000  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' ${params.Plink2_Options} --make-pgen 'erase-phase' --out \$name
   fi
  fi  
  """
}
//'@_#_\$r_\$a' is limited by Ref / Allele length although can be upped

