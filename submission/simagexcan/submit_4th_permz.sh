# ARGS1: config file middle name
# ARGS1: phenotype list

configname=$1
phenolist=$2

mkdir -p configs

for pheno in $(cat "${phenolist}"); do
  qsub -v CONFIGNAME="${configname}",GWASNAME="${pheno}" \
    -N ${configname}_$nametag run_simgx_4th_permz.qsub
done

