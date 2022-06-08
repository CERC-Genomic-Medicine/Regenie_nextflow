//Aim : Regenie parallelisation set-up

  input:
  tuple val(pheno_chunk_no), file(pheno_chunk), file(genotypes_file), file(sample_file), file(additional_file) from chunks_phenotypes.map { f -> [f.getBaseName().split('_')[1], f] }.combine(Channel.fromPath(params.genotypes_file).map(f -> f.getExtension() == "pgen" ? [f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), ""]))
  each file(covar_file) from Channel.fromPath(params.covar_file)


  output:
  tuple val(pheno_chunk_no), file(pheno_chunk), file("fit_bin${pheno_chunk_no}.master"), file("fit_bin${pheno_chunk_no}_*.snplist") into step1_l0_split mode flatten
  file "*.log" into step1_l0_logs

  publishDir "${params.OutDir}/step1_l0_logs", pattern: "*.log", mode: "copy"


  """
  if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
     echo "PGEN"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ "${params.CatCovar}" = "" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi

  if [ ${params.BinaryTrait} = "true" ]; then
     BinaryTrait="--bt"
  else
      BinaryTrait=""
      if [ ${params.RINT} = "true" ]; then
      BinaryTrait="--apply-rint"
      fi
  fi

  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \$BinaryTrait \
    --gz \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \$input \
    --out fit_bin_${pheno_chunk_no} \
    --split-l0 fit_bin${pheno_chunk_no},${params.njobs} \
    --threads ${params.Threads_S_10} \
    --lowmem
  """
}


process step_1_l1 {
  label "STEP_1_1"
  cache "lenient"
  scratch false

// Aim : Parallel Ridge Prediction

  input:
  tuple val(pheno_chunk_no), file(pheno_chunk), file(master), file(snplist), file(genotypes_file), file(sample_file), file(additional_file) from step1_l0_split.combine(Channel.fromPath(params.genotypes_file).map(f -> f.getExtension() == "pgen" ? [f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), ""]))
  each file(covar_file) from Channel.fromPath(params.covar_file)

  output:
  tuple val(pheno_chunk_no), file(pheno_chunk), file(master), file("*_l0_Y*") into step_1_l1
  file "*.log" into step1_l1_logs
 
  publishDir "${params.OutDir}/step1_l1_logs", pattern: "*.log", mode: "copy"

  """
  if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
     echo "PGEN"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ "${params.CatCovar}" = "" ]; then
     CovarCat=""
     echo "WORK"
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi

  if [ ${params.BinaryTrait} = "true" ]; then
     BinaryTrait="--bt"
  else 
      BinaryTrait=""
      if [ ${params.RINT} = "true" ]; then
      BinaryTrait="--apply-rint"
      fi
  fi

  i=${snplist.getSimpleName().split('_')[2].replaceFirst('^job', '')}
  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \$BinaryTrait \
    --gz \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --out fit_bin_${pheno_chunk_no}_\${i} \
    --run-l0 ${master},\${i} \
    --threads ${params.Threads_S_11} \
    --lowmem
  """
}


process step_1_l2 {
  label "STEP_1_2"
  cache "lenient"
  scratch false

// Aim : Gobal Prediction

  input:
  tuple val(pheno_chunk_no), file(pheno_chunk), file(master), file(predictions), file(genotypes_file), file(sample_file), file(additional_file) from step_1_l1.groupTuple(by: 0).map{ t -> [t[0], t[1][0], t[2][0], t[3].flatten()] }.combine(Channel.fromPath(params.genotypes_file).map(f -> f.getExtension() == "pgen" ? [f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), ""]))
  each file(covar_file) from Channel.fromPath(params.covar_file)

  output:       
  tuple val(pheno_chunk_no), file(pheno_chunk), file("fit_bin${pheno_chunk_no}_loco_pred.list"), file("*.loco.gz") into step1_l2_burden
  file "*.log" into step1_l2_logs

  publishDir "${params.OutDir}/step1_l2_logs", pattern: "*.log", mode: "copy"
  

  """
  if [ ${genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${genotypes_file.getBaseName()}"
  else
     input="--bgen ${genotypes_file} --sample ${sample_file}"
  fi

  if [ "${params.CatCovar}" = "" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi

  if [ ${params.BinaryTrait} ]; then
     BinaryTrait="--bt"
  else 
      BinaryTrait=""
      if [ ${params.RINT} = "true" ]; then
      BinaryTrait="--apply-rint"
      fi
  fi

  regenie \
    --step 1 \
    --loocv \
    --bsize ${params.Bsize} \$BinaryTrait \
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

// _________________________________________STEP2_________________________________________
process chunk_set {
   cache "lenient"
   scratch false
   executor "local"
   cpus 1

   input:
  file sets from Channel.fromPath(params.set_list)

   output:
   tuple val("${sets.getSimpleName()}"), file("${sets.getSimpleName()}_*.txt") into sets_chunks mode flatten

   """
  cut -f 1 ${sets} | split --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt -l ${params.SnpStep} - ${sets.getSimpleName()}_
   """
}
sets_chunks = sets_chunks.combine(Channel.fromPath(params.anno_file).map(f ->["${f.getSimpleName()}", f]), by: 0)
sets_chunks = sets_chunks.combine(Channel.fromPath(params.mask_def).map(f ->["${f.getSimpleName()}", f]), by: 0)
sets_chunks = sets_chunks.combine(Channel.fromPath(params.set_list).map(f ->["${f.getSimpleName()}", f]), by: 0)
sets_chunks = sets_chunks.combine(Channel.fromPath(params.gwas_genotypes_files).map(f -> f.getExtension() == "pgen" ? ["${f.getSimpleName()}", f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : ["${f.getSimpleName()}", f, file("${f.getParent()}/${f.getBaseName()}.sample"), f + ".bgi"]), by: 0)
sets_chunks.into { sets_chunks; sets_chunks_burden }
sets_chunks_burden.view()
// Decision Channel

Switch_non_burden_ch = (params.burden_testing & params.burden_testing_only ? Channel.empty():Channel.fromPath(params.gwas_genotypes_files))

process chunk_chromosomes {
   cache "lenient"
   scratch false
   executor "local"
   cpus 1

   input:
   file variants_file from Switch_non_burden_ch.map(f -> f.getExtension() == "pgen" ? file("${f.getParent()}/${f.getBaseName()}.pvar") : f + ".bgi")

   output:
   tuple val("${variants_file.getSimpleName()}"), file("${variants_file.getSimpleName()}_*.txt") into chromosome_chunks mode flatten

   """
   if [ ${variants_file.getExtension()} = "pvar" ]; then
      grep -v "^#" ${variants_file} | cut -f3
   else
      sqlite3 ${variants_file} "SELECT rsid FROM Variant"
   fi | split --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt -l ${params.SnpStep} - ${variants_file.getSimpleName()}_
   """
}


chromosome_chunks = chromosome_chunks.combine(Channel.fromPath(params.gwas_genotypes_files).map(f -> f.getExtension() == "pgen" ? ["${f.getSimpleName()}", f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : ["${f.getSimpleName()}", f, file("${f.getParent()}/${f.getBaseName()}.sample"), f + ".bgi"]), by: 0)


// _________________________________________STEP2_________________________________________


//___________________STEP 2 main ____________________________

process step_2_burden {
  label "Asscociation_testing_burden"
  cache "lenient"
  scratch false 

//Aim : Association testing

  input:
  tuple val(pheno_chunk_no), file(pheno_chunk), file(loco_pred_list), file(loco_pred), val(simple_name), file(set_chunk), file(Anno), file(mask), file(set), file(gwas_genotypes_file), file(samples_file), file(variants_file) from step1_l2_burden.combine(sets_chunks)
  each file(covar_file) from Channel.fromPath(params.covar_file)
  
  output:       
  file("*.regenie.gz") into summary_stats_burden mode flatten
  file "*.log" into step2_logs_burden
  
  publishDir "${params.OutDir}/step2_logs_burden", pattern: "*.log", mode: "copy"

  """
  if [ ${gwas_genotypes_file.getExtension()} = "pgen" ]; then
     input="--pgen ${gwas_genotypes_file.getBaseName()}"
  else
     input="--bgen ${gwas_genotypes_file} --sample ${samples_file}"
  fi

  if [ "${params.CatCovar}" = "" ]; then
     CovarCat=""
  else
      CovarCat="--catCovarList ${params.CatCovar}"
  fi
  if [ ${params.BinaryTrait} ]; then
    if [ ${params.firth} ]; then
      if [ ${params.approx} ]; then
        BinaryTrait='--bt --firth --approx --pThresh ${params.pThresh}'
      else
        BinaryTrait="--bt --firth --pThresh ${params.pThresh}"
      fi
    else
      BinaryTrait="--bt --SPA --pThresh ${params.pThresh}"
    fi
  else 
      BinaryTrait=""
      if [ ${params.RINT} = "true" ]; then
        BinaryTrait="--apply-rint"
      fi
  fi

  regenie \
    --step 2 \
    --gz \
    --loocv \
    --anno-file ${Anno} \
    --set-list ${set} \
    --mask-def ${mask} \
    --bsize ${params.Bsize} \$BinaryTrait \
    --phenoFile ${pheno_chunk} \
    --covarFile ${covar_file} \$CovarCat \
    \${input} \
    --test ${params.test} \
    --vc-tests ${params.test_Burden} \
    --out "${pheno_chunk_no}_${set_chunk.getSimpleName()}_assoc" \
    --pred ${loco_pred_list} \
    --extract-setlist ${set_chunk} \
    --threads ${params.Threads_S_2} \
    --lowmem
  """
}


//______________________MERGE______________________

process step_2_merge_burden {
  label "Merging_burden"
  cache "lenient"
  scratch false 

  //Aim : Natural Order Concatenated Association file (1 per phenotype)

  input:
  tuple val(pheno_name), file(summary) from summary_stats_burden.map{ t -> [t.baseName.split("_assoc_")[1], t] }.groupTuple()

  output:       
  file "*.txt.gz" into summary_stats_final_burden
  
  publishDir "${params.OutDir}/step2_result", pattern: "*.txt.gz", mode: "copy"

  """
  gzip -dc `find . -name "*.regenie.gz" -print -quit` | head -n1 | gzip -c > ${pheno_name}.txt.gz
  for f in `find . -name "*.regenie.gz" | sort -V`; do 
     gzip -dc \$f; 
  done | grep -v "^CHROM" | gzip -c >> ${pheno_name}.txt.gz
  """
}


