#!/usr/bin/env nextflow

/*
* AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2024
*/


process split_phenotypes {
   executor 'local'
   cache 'lenient'
  
   input:
   path phenotypes_file, stageAs: "phenotypes_file.txt"
  
   output:
   path "*.phenotype.txt"

   script:
   """
   # make sure phenotype file is tab-delimited
   cat ${phenotypes_file} | tr " " "\t" > temp_phenotypes_file.txt  

   n_cols=`head -n 1 temp_phenotypes_file.txt | wc -w` 
   
   for ((i=3;i<=\${n_cols};i++)); do
      name=`head -n 1 temp_phenotypes_file.txt | cut -f\${i}`
      cut -f 1,2,\${i} temp_phenotypes_file.txt > \${name}.phenotype.txt
   done
   """
}


process split_first_level_ridge_regression {
   executor 'local'
   cache 'lenient'

   container "${params.regenie_container}"

   input:
   tuple path(phenotypes_file), path(genotypes_file), path(sample_file), path(additional_file)
   each path(covariates_file)

   output:
   tuple path(phenotypes_file), path("${phenotypes_file.getName()}.master"), path("${phenotypes_file.getName()}_job*.snplist"), emit: out

   script:
   """
   if [ ${genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   if [ -z "${params.sex_specific}" ]; then
      sex_specific=""
   else
      sex_specific="--sex-specific ${params.sex_specific}"
   fi

   # This step doesn't need multi-threading because it just decides on distribution of the blocks between the machines for parallel processini
   regenie \
      --step 1 \
      --loocv \
      --skip-dosage-comp \
      --bsize ${params.block_size} \
      --gz \
      --phenoFile ${phenotypes_file} \
      --covarFile ${covariates_file} \${categorical_covariates} \
      \${input_genotypes} ${params.apply_rint} \${sex_specific} \
      --out ${phenotypes_file.getName()} \
      --split-l0 ${phenotypes_file.getName()},${params.n_ridge_regression_jobs} \
      --threads 1 \
      --lowmem
   """
}


process run_first_level_ridge_regression {
   cache 'lenient'

   cpus 8
   memory "16 GB"
   time "4h"

   container "${params.regenie_container}"

   input:
   tuple path(phenotypes_file), path(master_file), path(snplist_file), path(genotypes_file), path(sample_file), path(additional_file)
   each path(covariates_file)

   output:
   tuple val("${phenotypes_file.getName()}"), path("*_l0_Y*"), emit: out
   path "*.log", emit: log
 
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*_l0_Y*", mode: "copy"
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*.log", mode: "copy"
  
   script:
   """
   if [ ${genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   if [ -z "${params.sex_specific}" ]; then
      sex_specific=""
   else
      sex_specific="--sex-specific ${params.sex_specific}"
   fi

   i=`echo "${snplist_file}" | sed -nr 's/.*_job([1-9][0-9]*).snplist/\\1/p'`
   regenie \
      --step 1 \
      --loocv \
      --skip-dosage-comp \
      --bsize ${params.block_size} \
      --gz \
      --phenoFile ${phenotypes_file} \
      --covarFile ${covariates_file} \${categorical_covariates} \
      \${input_genotypes} ${params.apply_rint} \${sex_specific} \
      --out ${phenotypes_file.getName()}_job\${i}_first_level_ridge_regression \
      --run-l0 ${master_file},\${i} \
      --threads 8 \
      --lowmem
   """
}


process run_second_level_ridge_regression {
   cache 'lenient'

   cpus 8
   memory "16GB"
   time "1h"

   container "${params.regenie_container}"

   input:
   tuple path(phenotypes_file), path(master_file), path(first_level_predictions), path(genotypes_file), path(sample_file), path(additional_file)
   each path(covariates_file)

   output:
   tuple val("${phenotypes_file.getName()}"), path("*_pred.list"), path("*.loco.gz"), emit: out
   path "*.log", emit: logs

   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*_pred.list", mode: "copy"
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*.loco.gz", mode: "copy"
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*.log", mode: "copy"

   script:
   """
   if [ ${genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   if [ -z "${params.sex_specific}" ]; then
      sex_specific=""
   else
      sex_specific="--sex-specific ${params.sex_specific}"
   fi

   regenie \
      --step 1 \
      --loocv \
      --skip-dosage-comp \
      --bsize ${params.block_size} \
      --gz \
      --phenoFile ${phenotypes_file} \
      --covarFile ${covariates_file} \${categorical_covariates} \
      \${input_genotypes} ${params.apply_rint}  \${sex_specific} \
      --out ${phenotypes_file.getName()} \
      --run-l1 ${master_file} \
      --keep-l0 \
      --threads 8 \
      --use-relative-path \
      --lowmem
   """
}


process run_all_ridge_regressions {
   cache 'lenient'

   errorStrategy {
      def delay = Math.pow(2, task.attempt) * 300 as long
      println "Retrying ${task.process} attempt ${task.attempt} after ${delay}ms"
      sleep(delay)
      return 'retry'
   }

   cpus 8
   memory "16GB"
   time "1h"

   container "${params.regenie_container}"

   input:
   tuple path(phenotypes_file), path(genotypes_file), path(sample_file), path(additional_file)
   each path(covariates_file)

   output:
   tuple val("${phenotypes_file.getName()}"), path("*_pred.list"), path("*.loco.gz"), emit: out
   path "*.log", emit: logs

   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*_pred.list", mode: "copy"
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*.loco.gz", mode: "copy"
   publishDir "${params.output_dir}/genomic_predictions/${phenotypes_file.getName()}/", pattern: "*.log", mode: "copy"
 
   script:
   """
   if [ ${genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   if [ -z "${params.sex_specific}" ]; then
      sex_specific=""
   else
      sex_specific="--sex-specific ${params.sex_specific}"
   fi

   # Use `--strict` option when analyzing one phenotype at a time to avoid imputation of missing phenotype values (i.e. keep only individuals with non-missing phenotypes).
   regenie \
      --step 1 \
      --loocv \
      --skip-dosage-comp \
      --bsize ${params.block_size} \
      --gz \
      --phenoFile ${phenotypes_file} --strict \
      --covarFile ${covariates_file} \${categorical_covariates} \
      \${input_genotypes} ${params.apply_rint} \${sex_specific} \
      --out ${phenotypes_file.getName()} \
      --threads 8 \
      --lowmem \
      --use-relative-path
   """
}


process chunk_chromosomes {
   cache "lenient"
   executor "local"

   input:
   path variants_file

   output:
   tuple val("${variants_file.getBaseName()}"), path("${variants_file.getBaseName()}_*.txt")

   """
   # Get the list of variant IDs
   if [ ${variants_file.getExtension()} = "pvar" ]; then
         grep -v "^#" ${variants_file} | cut -f3
   else
         sqlite3 ${variants_file} "SELECT rsid FROM Variant"
   fi | gzip -c > sites.txt.gz

   # Choose the maximal number of equal chunks such that they don't exceed the maximal number of variants per chunk
   n_chunks=1
   n_sites=`gzip -dc sites.txt.gz | split -n l/1/\${n_chunks} | wc -l`
   while [ \${n_sites} -gt ${params.gwas_chunk_size} ]; do
         n_chunks=\$((n_chunks+1))
         n_sites=`gzip -dc sites.txt.gz | split -n l/1/\${n_chunks} | wc -l`
   done

   # Final chunking
   gzip -dc sites.txt.gz | split --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt -n l/\${n_chunks} - ${variants_file.getBaseName()}_
   """
}


process run_association_tesing {
   cache "lenient"
   //scratch false

   maxRetries 3
   errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }

   cpus 6
   // memory "24GB"
   memory { 24.GB + ((task.attempt ?: 1) - 1) * 8.GB }
   time "1h"

   container "${params.regenie_container}"

   input:
   tuple path(phenotypes_file), path(loco_pred_list), path(loco_pred), path(chromosome_chunk), path(gwas_genotypes_file), path(sample_file), path(variants_file)
   each path(covariates_file)
  
   output:       
   tuple val("${chromosome_chunk.getBaseName()}"), path("*.regenie.gz"), emit: out
   path("*.log"), emit: logs
  
   publishDir "${params.output_dir}/GWAS/logs", pattern: "*.log", mode: "copy"

   """
   if [ ${gwas_genotypes_file.getExtension()} = "pgen" ]; then
      input_genotypes="--pgen ${gwas_genotypes_file.getBaseName()}"
   else
      input_genotypes="--bgen ${gwas_genotypes_file} --sample ${sample_file}"
   fi

   if [ -z "${params.categorical_covariates}" ]; then
      categorical_covariates=""
   else
      categorical_covariates="--catCovarList ${params.categorical_covariates}"
   fi

   if [ ${params.split_phenotypes} = true ]; then
      strict="--strict" # Exclude individuals with any missing phenotypes.
   else
      strict=""
   fi

   if [ -z "${params.gene_by_sex_interaction}" ]; then
      interaction_variable=""
   else
      interaction_variable="--interaction ${params.gene_by_sex_interaction}"
   fi

   if [ -z "${params.sex_specific}" ]; then
      sex_specific=""
   else
      sex_specific="--sex-specific ${params.sex_specific}"
   fi

   regenie \
    --step 2 \
    --gz \
    --loocv \
    --skip-dosage-comp \
    --bsize ${params.block_size} \
    --phenoFile ${phenotypes_file} \${strict} \
    --covarFile ${covariates_file} \${categorical_covariates} \${interaction_variable} \
    \${input_genotypes} ${params.apply_rint} \${sex_specific} \
    --out "${chromosome_chunk.getBaseName()}" \
    --pred ${loco_pred_list} \
    --extract ${chromosome_chunk} \
    --threads 8 \
    --lowmem
   """
}


process merge_association_results {
   cache "lenient"
   //scratch false

   cpus 1
   memory "2GB"
   time "1h"

   input:
   tuple val(filename), path(chunked_association_results)

   output:       
   path "${filename}", emit: out
  
   publishDir "${params.output_dir}/GWAS/results", pattern: "${filename}", mode: "copy"

   """
   gzip -dc `find . -name "*_${filename}" -print -quit` | head -n1 | gzip -c > ${filename}
   for f in `find . -name "*_${filename}" | sort -V`; do
      gzip -dc \$f; 
   done | grep -v "^CHROM" | gzip -c >> ${filename}
   """
}

workflow {
   if (params.sex_specific && params.gene_by_sex_interaction) {
      println "Error: male-only or female-only analyses can't be run together with the gene-by-sex interaction."
      return 1
   }

   // Load phenotypes
   if (params.split_phenotypes == true) {
      phenotypes = split_phenotypes(Channel.fromPath(params.phenotypes_file)).flatten()
   } else {
      phenotypes = Channel.fromPath(params.phenotypes_file)
   }
   //phenotypes.view()

   // Load covariates
   covariates = Channel.fromPath(params.covariates_file)
   //covariates.view()

   // Genomic predictions
   if (params.genomic_predictions_files) {
      // If specified, then load pre-computed genomic predictions for phenotypes
      genomic_predictions = Channel.fromPath(params.genomic_predictions_files).map(it -> [it.getName().replaceAll(/_pred.list$/, ""), it, files(it.toString().replaceAll(/_pred.list$/, "_*.loco.gz"), checkIfExists: true)])
   } else {
      // Load pruned genotypes
      pruned_genotypes = Channel.fromPath(params.pruned_genotypes_file).map(f -> f.getExtension() == "pgen" ? [f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : [f, file("${f.getParent()}/${f.getBaseName()}.sample"), ""])
      //pruned_genotypes.view()

      // If multiple phenotypes are split and analyzed one by one, then each phenotype is analyzed in a single compute job/machine.
      // If multiple phenotypes are analyzed in bulk, then compute is distributed to multiple jobs/machines by chunking chromosome and the number of jobs is controlled by `params.n_ridge_regression_jobs`.
      if (params.split_phenotypes == true) {
         genomic_predictions = run_all_ridge_regressions(phenotypes.combine(pruned_genotypes), covariates).out
      } else {
         first_level_ridge_regression_jobs = split_first_level_ridge_regression(phenotypes.combine(pruned_genotypes), covariates).out
         //first_level_ridge_regression_jobs.view()

         scattered_first_level_predictions = run_first_level_ridge_regression(first_level_ridge_regression_jobs.transpose().combine(pruned_genotypes), covariates).out
         //scattered_first_level_predictions.view()
     
         gathered_first_level_predictions = first_level_ridge_regression_jobs.map(it -> [it[0].getName(), it[0], it[1]]).join(
            scattered_first_level_predictions.groupTuple(by: [0]), by: [0], failOnMismatch: true, failOnDuplicate: true).map(it -> [it[1], it[2], it[3].flatten()])
         //gathered_first_level_predictions.view()

         genomic_predictions = run_second_level_ridge_regression(gathered_first_level_predictions.combine(pruned_genotypes), covariates).out
      }
   }
   //genomic_predictions.view()
 
   // Load GWAS genotypes:
   gwas_genotypes = Channel.fromPath(params.gwas_genotypes_files).map(f -> f.getExtension() == "pgen" ? ["${f.getBaseName()}", f, file("${f.getParent()}/${f.getBaseName()}.psam"), file("${f.getParent()}/${f.getBaseName()}.pvar")] : ["${f.getBaseName()}", f, file("${f.getParent()}/${f.getBaseName()}.sample"), f + ".bgi"])
   //gwas_genotypes.view()

   // Split gwas variants into chunks for parallel processing based on pvar (PLINK) or bgi (BGEN) files. 
   gwas_chunks = chunk_chromosomes(gwas_genotypes.map(it -> it[3])).transpose().combine(gwas_genotypes, by: [0]).map(it -> it.drop(1))
   //gwas_chunks.view()
  
   // Run GWAS
   association_results_by_chunk = run_association_tesing(
      phenotypes.map(it -> [it.getName(), it]).join(genomic_predictions, by: [0], failOnMismatch: true, failOnDuplicate: true).map(it -> it.drop(1)).combine(gwas_chunks), 
      covariates).out

   merge_association_results(association_results_by_chunk.transpose().map(it -> [it[1].getName().replaceAll(/${it[0]}_/, ""), it[1]]).groupTuple(by: [0]))
}
