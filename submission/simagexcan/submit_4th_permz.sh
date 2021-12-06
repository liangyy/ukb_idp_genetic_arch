# ARGS1: config file middle name
# ARGS1: phenotype list

configname=$1
phenolist=$2

mkdir -p configs

for pheno in $(cat "${phenolist}"); do
  if [[ $pheno != "Neuroticism_CTG" && $pheno != "AlcDep_PGC_2018" ]]; then
    extra="effect_size=effect_size"
  else
    extra=""
  fi
  ff="logs/run_simgx_4th_permz_${configname}_${pheno}.log"
  if [[ -f "${ff}" ]]; then
    tmp=$(cat "${ff}" | tail -n 2 | grep "failed\|kill\|Errno\|File" | wc -l)
    if [[ $tmp != 0 ]]; then
      qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}",EXTRA="${extra}" \
        -N ${configname}_${pheno} run_simgx_4th_permz.qsub
    fi
  else
    qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}",EXTRA="${extra}" \
      -N ${configname}_${pheno} run_simgx_4th_permz.qsub
  fi
  
done

