# ARGS1: nametag
# ARGS2: phenotype_list (config middle name)
# ARGS3: processor per node
nametag=$1
pheno_list=$2
ppn=$3

mkdir -p configs

cat config.$pheno_list.yaml | sed "s#PLACEHOLDER#$nametag#g" > configs/config.${pheno_list}_$nametag.yaml

qsub -v NAME=${pheno_list}_$nametag,PPN=$ppn -l nodes=1:ppn=$ppn -N ${pheno_list}_$nametag run_simgx_3rd.qsub
