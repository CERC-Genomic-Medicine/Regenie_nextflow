
process chunk_phenotype {
  label "chunk"
  executor "local"
  
  publishDir params.OutDir
  input :
// channel from path for phename
  output:
  file "chunk_*_phe.txt" into chunks_phenotypes,chunks_phenotypes_l0, chunks_phenotypes_l1, chunks_phenotypes_l2,numero  mode flatten
  stdout into Q
  """
sed -i.bak \$'s/\t/ /g' ${params.InDir}"/"${params.PheName}
Nb_PHENO=\$((\$(head ${params.InDir}"/"${params.PheName} -n 1 | wc -w ) - 2))
val=\$((\$Nb_PHENO/${params.PheStep}))
if [ \$val > 1 ]
then
for ((Q=1;Q<=\$val;Q++)); do
cat ${params.InDir}/${params.PheName} | cut -f 1,2,\$((( \$Q - 1) * ${params.PheStep} + 3 ))-\$(((\$Q * ${params.PheStep}) + 2)) -d " " > chunk_"\$Q"_phe.txt
done
fi
if [ \$val <=  1 ]
then
cat ${params.InDir}/${params.PheName}  | cut -f 1,2,\$((\$val * ${params.PheStep} + 3))-\$((\$Nb_PHENO + 2)) -d " " > chunk_"\$Q"_phe.txt
fi
echo \$Q
	"""
}

(Phenos , Phenos_l0, Phenos_l1, Phenos_l2, Phenos_s2) =numero.flatMap {n ->  n.getBaseName().split('_')[1] }.into(4)





process step1_l0 {
	label "STEP_1_0"
	containerOptions "-B ${params.InDir}:$HOME/input"
	cpus = 2
  time = "12h"
  memory = "6GB"

	input:
file(chunk) from chunks_phenotypes 
val(pheno) from Phenos
	output:
 file '*.master' into masterFiles,masterFiles_l1,masterFiles_l2
 file '*.log' into logs
 file "*.snplist" into snplists,snplists_l1,snplists_l2
"""
	regenie \
    --step 1 \
    --loocv \
    --phenoFile ${chunk} \
    --bsize ${params.Bsize} \
    --gz \
    --bgen \$HOME/input/${params.bfile} \
    --out test_bin_{$pheno} \
    --split-l0 fit_bin${pheno},${params.njobs} \
    --threads 1 \
    --extract \$HOME/input/qc_pass.snplist  \
    --force-step1 ${params.options}
	"""
}


def range = 1..params.njobs
iter = Channel.from(range.by(1))

(iter,iter_s2)=iter.into(2)

process step_1_l1 {
	label "STEP_1_1"
	containerOptions "-B ${params.InDir}:$HOME/input"
  cpu=1
	input:
    file(master) from masterFiles
     each i from iter 
      file(chunk) from chunks_phenotypes_l0
      val(pheno) from Phenos_l0
      file(snplist) from snplists.collect()
      file(log) from logs
 


	output:
 file("fit_bin*") into splited_files,splited_files_l2
 file("*.done") into omega
 
"""
	regenie \
    --step 1 \
    --loocv \
    --phenoFile ${chunk} \
    --bsize ${params.Bsize} \
    --gz \
    --bgen \$HOME/input/${params.bfile} \
    --out test_bin_${pheno} \
    --run-l0 ${master},${i} \
    --extract \$HOME/input/qc_pass.snplist  \
    --threads 1 ${params.options}
    
    touch ${i}.done
	"""
}



process step_1_l2 {
	label "STEP_1_2"
	containerOptions "-B ${params.InDir}:$HOME/input"
 	cpus 1
  
  publishDir params.OutDir

    input:
  file m from masterFiles_l1
  file(snp) from splited_files.collect() //doesn't wait for all iteration the last one to be finished
  val(pheno) from Phenos_l1
  file(chunk) from chunks_phenotypes_l1

    output:       
  set val(pheno), file(chunk),  file("test_bin_*_pred.list") into pred_list
  file "*.loco.gz" into preds,pred_names 
     """
	regenie \
    --step 1 \
    --loocv \
    --covarFile \$HOME/input/${params.CovarName} \
    --phenoFile ${chunk} \
    --bsize ${params.Bsize} \
    --gz \
    --bgen $HOME/input/${params.bfile} \
    --out test_bin_${pheno} \
    --run-l1 ${m} \
    --keep-l0 \
    --threads 1 \
    --extract \$HOME/input/qc_pass.snplist  \
    --force-step1 ${params.options}
	"""
}



process step_2 {
	label "STEP_2"
	containerOptions "-B ${params.InDir}:$HOME/input"
 	cpus 1
    input:
  set val(pheno), file(chunk), file(pred_l) from pred_list
  each i from iter_s2
  file(snp) from snplists_l1.collect()
  
  file preds from preds.collect()

   output:       
  file "*.regenie" into regenies_split, regenies_split_names mode flatten
     """
    regenie \
    --step 2 \
    --phenoFile ${chunk} \
    --bsize ${params.Bsize} \
    --bgen $HOME/input/${params.bfile} \
    --out step2_${i}.${pheno}.r\
    --pred ${pred_l} \
    --threads 1 \
    --extract ${snp} ${params.options}
	"""
}

pheno_names = regenies_split_names.flatMap(n-> n.getBaseName().split('.r')[1]).unique()



process concat {
  executor "local"
	label "Concat"
 
 publishDir params.OutDir

	input :
	val(ph) from pheno_names
	file(regenie) from regenies_split.collect()

	output :
	file("*.regenie.gz") into ENDING

	"""
	cat step2_1.*.r${ph}.regenie > regenie.out
	for ((J=2; J<=${params.njobs}; J++))
		do
		sed 1d step2_"\$J".*.r${ph}.regenie >> regenie.out
		done
	zip ${ph}.regenie.gz regenie.out
	"""
	}
