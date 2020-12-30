# ARGS1: nametag
# ARGS2: phenotype_list (config middle name)
# ARGS3: processor per node
nametag=$1
pheno_list=$2
ppn=$3
cat config.$pheno_list.yaml | sed "s#PLACEHOLDER##g" > configs/config.${pheno_list}_$nametag.yaml

qsub -v PPN=$ppn,NAME=${pheno_list}_$nametag -N ${pheno_list}_$nametag run_simgx_2nd.qsub
