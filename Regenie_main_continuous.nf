process CHUNK_PHENOTYPE {
  label 'chunk'
  executor 'local'
  cache 'lenient'
  
  input:
  path pheno_file
  
  output:
  path "chunk_*_phe.txt"

  publishDir "${params.OutDir}/chunked_pheno", pattern: "chunk_*_phe.txt", mode: "copy"

  script:
  """
    # make sure phenotype file is tab-delimited
    cat ${pheno_file} | tr " " "\t" > temp_pheno_file.txt  

    Nb_PHENO=\$((\$(head -n 1 temp_pheno_file.txt | wc -w ) - 2)) 
    val=\$((\$Nb_PHENO/${params.PheStep}))
    mod=\$((\$Nb_PHENO%${params.PheStep}))
    if [[ \$val > 0 ]]; then
        for ((Q=1;Q<=\$val;Q++)); do
            cut -f 1,2,\$((( \$Q - 1) * ${params.PheStep} + 3 ))-\$(((\$Q * ${params.PheStep}) + 2)) temp_pheno_file.txt > chunk_\${Q}_phe.txt
        done
        if [[ \$mod != 0 ]]; then
            cut -f 1,2,\$((( \$Q - 1) * ${params.PheStep} + 3 ))-\$(\$Nb_PHENO + 3) temp_pheno_file.txt > chunk_\${Q}_phe.txt
        fi
    else
        cp temp_pheno_file.txt chunk_1_phe.txt
    fi
  """
}

process STEP1_L0 {
  label 'STEP_1_0'
  cache 'lenient'
  scratch true

  input:
  tuple val(pheno_chunk_no), path(pheno_chunk), path(genotypes_file), path(sample_file), path(additional_file)
  each path(covar_file)

  output:
  tuple val(pheno_chunk_no), path(pheno_chunk), path("fit_bin${pheno_chunk_no}.master"), path("fit_bin${pheno_chunk_no}_*.snplist"), emit: step1_l0
  path "*.log", emit: step1_l0_logs

  publishDir "${params.OutDir}/step1_l0/step1_l0_logs", pattern: "*.log", mode: "copy"
  publishDir "${params.OutDir}/step1_l0/step1_l0_${pheno_chunk_no}", pattern: "fit_bin${pheno_chunk_no}_*.snplist", mode: "copy"

  script:
  """
if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ -z "${params.CatCovar}" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi

  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \
    --gz \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --out fit_bin_${pheno_chunk_no} \
    --split-l0 fit_bin${pheno_chunk_no},${params.njobs} \
    --threads ${params.Threads_S_10} \
    --lowmem
  """
}

process STEP_1_L1 {
  label 'STEP_1_1'
  cache 'lenient'
  scratch false

  input:
  tuple val(pheno_chunk_no), path(pheno_chunk), path(master), path(snplist),val(run), path(genotypes_file), path(sample_file), path(additional_file)
  each path(covar_file)

  output:
  tuple val(pheno_chunk_no), path(pheno_chunk), path(master),path("*_l0_Y*"), emit: step_1_l1_out
  path "*.log", emit: step1_l1_logs

  publishDir "${params.OutDir}/step1_l1/step1_l1_logs", pattern: "*.log", mode: "copy"
  publishDir "${params.OutDir}/step1_l1/step1_l1_chunk_${pheno_chunk_no}/", pattern: "*_l0_Y*", mode: "copy"

  script:
  """
  if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ -z "${params.CatCovar}" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi


  i=${run}
  echo \$i
  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \
    --gz \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --out \${i} \
    --run-l0 ${master},\${i} \
    --threads ${params.Threads_S_11} \
    --lowmem
  """
}

process STEP_1_L2 {
  label 'STEP_1_2'
  cache 'lenient'
  scratch false

  input:
  tuple val(pheno_chunk_no), path(pheno_chunk), path(master),path(predictions), path(genotypes_file), path(sample_file), path(additional_file)
  each path(covar_file)

  output:
  tuple val(pheno_chunk_no), path(pheno_chunk), path("fit_bin${pheno_chunk_no}_loco_pred.list"), path("*.loco.gz"), emit: step1_l2_out
  path "*.log", emit: step1_l2_logs

  publishDir "${params.OutDir}/step1_l2/step1_l2_chunk_${pheno_chunk_no}", pattern: "*.loco.gz", mode: "copy"
  publishDir "${params.OutDir}/step1_l2/step1_l2_chunk_${pheno_chunk_no}", pattern: "fit_bin${pheno_chunk_no}_loco_pred.list", mode: "copy"
  publishDir "${params.OutDir}/step1_l2/logs", pattern: "*.log", mode: "copy"
  script:
  """
  if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ -z "${params.CatCovar}" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi

  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \
    --gz \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --out fit_bin${pheno_chunk_no}_loco \
    --run-l1 ${master} \
    --keep-l0 \
    --threads ${params.Threads_S_12} \
    --use-relative-path \
    --lowmem
    """
}

process chunk_chromosomes {
   cache "lenient"
   scratch false
   executor "local"
   cpus 1
   label 'SNP_chunk'

   input:
   path(variants_file)

   output:
   path("${variants_file.getSimpleName()}_*.txt"), emit: chromosome_chunks

   publishDir "${params.OutDir}/step2/Variant_Chunk", pattern: "*.txt", mode: "copy"

   """
   if [ ${variants_file.getExtension()} = "pvar" ]; then
      grep -v "^#" ${variants_file} | cut -f3
   else
      sqlite3 ${variants_file} "SELECT rsid FROM Variant"
   fi | split --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt -l ${params.SnpStep} - ${variants_file.getSimpleName()}_
   """
}


//___________________STEP 2 main ____________________________

process step_2 {
  label "Asscociation_testing"
  cache "lenient"
  scratch false


//Aim : Association testing

  input:
  tuple val(pheno_chunk_no), path(pheno_chunk), path(loco_pred_list), path(loco_pred), val(simple_name), path(chromosome_chunk), path(gwas_genotypes_file), path(samples_file), path(variants_file)
  each path(covar_file)
  
  output:       
  path("*.regenie.gz"), emit: summary_stats
  path("*.log"), emit: step2_logs
  path("*.list"), optional: true
  

  publishDir "${params.OutDir}/step2/result/", pattern: "*.regenie.gz", mode: "copy"
  publishDir "${params.OutDir}/step2/logs", pattern: "*.log", mode: "copy"

  """
  if [ ${gwas_genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${gwas_genotypes_file.getBaseName()}"
  else
     input="--bgen ${gwas_genotypes_file} --sample ${samples_file}"
  fi

  if [ -z "${params.CatCovar}" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi


  regenie \
    --step 2 \
    --gz \
    --loocv \
    --bsize ${params.Bsize} \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --out "${pheno_chunk_no}_${chromosome_chunk.getSimpleName()}_assoc" \
    --pred ${loco_pred_list} \
    --extract ${chromosome_chunk} \
    --threads ${params.Threads_S_2} \
    --lowmem
    """
}


//______________________MERGE______________________
process step_2_merge {
  label "Merging"
  cache "lenient"
  scratch false

  //Aim : Natural Order Concatenated Association file (1 per phenotype)

  input:
  tuple val(pheno_name), path(summary)

  output:       
  path "*.txt.gz", emit: summary_stats_final
  

  publishDir "${params.OutDir}/step2/summary", pattern: "*.txt.gz", mode: "copy"

  """
  gzip -dc `find . -name "*.regenie.gz" -print -quit` | head -n1 | gzip -c > ${pheno_name}.txt.gz
  for f in `find . -name "*.regenie.gz" | sort -V`; do
     gzip -dc \$f; 
  done | grep -v "^CHROM" | gzip -c >> ${pheno_name}.txt.gz
  """
}

workflow {
  //Input files
     Common_LD_pruned_variant = Channel.fromPath(params.genotypes_file).map(f -> f.getExtension() == "pgen" ? [f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), ""])
     GWAS_variant = Channel.fromPath(params.gwas_genotypes_files).map(f -> f.getExtension() == "pgen" ? ["${f.getSimpleName()}".split('_')[0], f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : ["${f.getSimpleName()}".split('_')[1], f, file("${f.getParent()}/${f.getBaseName()}.sample"), f + ".bgi"])
     Covariant=Channel.fromPath(params.covar_file)

  //chunking Phenotype
     chunks = CHUNK_PHENOTYPE(Channel.fromPath(params.pheno_file))
  //chunking SNP
     SNP_chunks = chunk_chromosomes(Channel.fromPath(params.gwas_genotypes_files).map(f -> f.getExtension() == "pgen" ? file("${f.getParent()}/${f.getBaseName()}.pvar") : f + ".bgi"))
     S=SNP_chunks.chromosome_chunks.flatten()
     T=S.map(t -> ["${t.getSimpleName()}".split('_')[0],t])
     SNP_chunk = T.combine(GWAS_variant, by: 0)

  //Regenie Step 1
     step1_L0_input = chunks.flatten().map { f -> [f.getBaseName().split('_')[1], f] }

     S1_L0=STEP1_L0(step1_L0_input.combine(Common_LD_pruned_variant), Covariant )
  
  // Divides each S1_L0 SNP list into a process
  // Includes the modeling data
     jobs_S1=Channel.from( 1..params.njobs )
     S1_L1_input=S1_L0.step1_l0.combine(jobs_S1)
     S1_L1=STEP_1_L1(S1_L1_input.combine(Common_LD_pruned_variant), Covariant)

  //Regroup the output of S1_L1 per phenotype_chuck
     STEP_1_L1_reformated = S1_L1.step_1_l1_out.groupTuple(by: 0).map{ t -> [t[0], t[1][0], t[2][0], t[3].flatten()] }

     S1 = STEP_1_L2(STEP_1_L1_reformated.combine(Common_LD_pruned_variant), Covariant)

  // Regenie Step 2
     S2_input = S1.step1_l2_out.combine(SNP_chunk)

    S2 = step_2(S2_input,Covariant)
    S2.summary_stats.flatten().map{ t -> [t.baseName.split('.')] }.view()
    S2_groups = S2.summary_stats.flatten().map{ t -> [t.baseName.split('_')[5],t] }.groupTuple()
  // Step 2 Merged
     S2_Merged = step_2_merge(S2_groups)

}
