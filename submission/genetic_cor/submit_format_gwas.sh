# ARGS1: pheno list
# ARGS2: gwas name tag

for pheno in `cat $1`
do
  qsub -v NAME=$2,GWASNAME=$pheno -N ${2}_$pheno format_gwas.qsub
done
