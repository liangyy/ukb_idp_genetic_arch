# ARGS1: pheno list
# ARGS2: gwas name tag

for pheno in `cat $1`
do
  ff=logs/format_gwas_${2}_$pheno.log
  if [[ -f $ff ]]
  then
    e=`cat $ff | tail -n 1 | grep 'failed\|kill\|Errno\|found' | wc -l` 
    if [[ $e = 1 ]]
    then
      qsub -v NAME=$2,GWASNAME=$pheno -N ${2}_$pheno format_gwas.qsub
    fi
  else
    qsub -v NAME=$2,GWASNAME=$pheno -N ${2}_$pheno format_gwas.qsub
  fi
done
