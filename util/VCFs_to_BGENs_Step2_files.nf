//_____________________Description : parallel conversion and filtering of VCF file to BGEN file

//unziping of exclusion region (ie low complexity region and or repeat sequence)
process Unziping {
  label "Unziping BED"
  cache "lenient"
  scratch true
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

//Conversion of each VCF/BCF to BGEN 
process Plinked {
  label "BGEN_generation"
  cache "lenient"
  scratch true
    
  input :
  tuple file(vcf), file(vcf_index) from Channel.fromPath(params.VCF_files).map{ vcf -> [ vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".tbi") ] } // accepts VCF or BCF
  each file(Exclude) from excluded_unzip

  output:
  file "*.bgen" into BGEN_file mode flatten
      publishDir "${params.OutDir}/BGEN", pattern: "*.bgen", mode: "copy"
  file "*.sample" into sample_file mode flatten
      publishDir "${params.OutDir}/BGEN", pattern: "*.sample", mode: "copy"
  file "*.bgi" into BGI_file mode flatten
      publishDir "${params.OutDir}/BGEN", pattern: "*.bgi", mode: "copy"
    
  """
  if [ ${vcf.getExtension()} = "bcf"  ];then
    name=${vcf.getName().replaceAll('.bcf$', '')}
   if [ ${params.dosageFD} = "none"  ];then
    plink2 --bcf ${vcf}  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' --new-id-max-allele-len 1000 --exclude ${Exclude} --export vcf-4.2 ref-first --out \$name
    bash ${workflow.scriptFile.getParent()}/Scripts/dephasing.sh -f "\$name".vcf
    qctool -g "\$name".vcf -filetype vcf  -bgen-bits 8 -og "\$name".bgen -os "\$name".sample
    bgenix -index -g "\$name".bgen 
   else
    plink2 --bcf ${vcf} 'dosage=${params.dosageFD}'  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' --new-id-max-allele-len 1000 --exclude ${Exclude} --export vcf-4.2 ref-first 'vcf-dosage='${params.dosageFD} --out \$name
    bash ${workflow.scriptFile.getParent()}/Scripts/dephasing.sh -f "\$name".vcf
    qctool -g "\$name".vcf -filetype vcf  -vcf-genotype-field ${params.dosageFD} -bgen-bits 8 -og "\$name".bgen -os "\$name".sample
    bgenix -index -g "\$name".bgen    
   fi
  else
   name=${vcf.getName().replaceAll('.vcf.gz$', '')}
   if [ ${params.dosageFD} = "none"  ];then
    plink2 --vcf ${vcf}  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' --new-id-max-allele-len 1000  --exclude ${Exclude} --export vcf-4.2 ref-first --out \$name
    bash ${workflow.scriptFile.getParent()}/Scripts/dephasing.sh -f "\$name".vcf
    qctool -g "\$name".vcf -filetype vcf  -bgen-bits 8 -og "\$name".bgen -os "\$name".sample
    bgenix -index -g "\$name".bgen 
   else
    plink2 --vcf ${vcf} 'dosage=${params.dosageFD}'  --max-alleles 2 --set-all-var-ids '@_#_\$r_\$a' --new-id-max-allele-len 1000 --exclude ${Exclude} --export vcf-4.2 ref-first 'vcf-dosage='${params.dosageFD} --out \$name
    bash ${workflow.scriptFile.getParent()}/Scripts/dephasing.sh -f "\$name".vcf
    qctool -g "\$name".vcf -filetype vcf  -vcf-genotype-field ${params.dosageFD} -bgen-bits 8 -og "\$name".bgen -os "\$name".sample
    bgenix -index -g "\$name".bgen 
   fi
  fi  
  """
}
//'@_#_\$r_\$a' is limited by Ref / Allele length although can be up
