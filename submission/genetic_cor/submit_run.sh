# ARGS1: pheno list
# ARGS2: gwas name tag
# ARGS3: IDP tag
# ARGS4: IDP list

for pheno in `cat $1`
do
  ff=logs/run_${2}_${pheno}_${3}.log
  if [[ -f $ff ]]
  then
    e=`cat $ff | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l` 
    if [[ $e = 1 ]]
    then
      qsub -v NAME=$2,GWASNAME=$pheno,IDP=$3,IDPLIST=$4 -N ${2}_$pheno_${3} run.qsub
    fi
  else
    qsub -v NAME=$2,GWASNAME=$pheno,IDP=$3,IDPLIST=$4 -N ${2}_$pheno_${3} run.qsub
  fi
done
