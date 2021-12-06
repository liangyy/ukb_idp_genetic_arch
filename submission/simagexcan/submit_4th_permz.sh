# ARGS1: config file middle name
# ARGS1: phenotype list

configname=$1
phenolist=$2

mkdir -p configs

for pheno in $(cat "${phenolist}"); do
  ff="logs/run_simgx_4th_permz_${configname}_${pheno}.log"
  if [[ -f "${ff}" ]]; then
    tmp=$(cat "${ff}" | tail -n 1 | grep "failed\|kill\|Errno\|File" | wc -l)
    if [[ $tmp = 1 ]]; then
      echo qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}" \
        -N ${configname}_${pheno} run_simgx_4th_permz.qsub
    fi
  else
    echo qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}" \
      -N ${configname}_${pheno} run_simgx_4th_permz.qsub
  fi
done

