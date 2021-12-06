# ARGS1: config file middle name
# ARGS1: phenotype list

configname=$1
phenolist=$2

mkdir -p configs

for pheno in $(cat "${phenolist}"); do
  if [[ $pheno == "Neuroticism_CTG" || $pheno == "AlcDep_PGC_2018" ]]; then
    extra="zscore=zscore"
  fi
  echo qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}",EXTRA="" \
    -N ${configname}_$nametag run_simgx_4th_permz.qsub
done

