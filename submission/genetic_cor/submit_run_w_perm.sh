# ARGS1: pheno id
# ARGS2: gwas name tag
# ARGS3: IDP tag
# ARGS4: IDP list
# ARGS5: nrepeat 
# ARGS6: seed
# ARGS7: LD block

pheno=$1
nrepeat=$5

# run w/o permutation
idpname="${3}"
ff=logs/run_${2}_${pheno}_${idpname}.log
if [[ -f $ff ]]
then
  e=`cat $ff | tail -n 1 | grep 'failed\|kill\|Errno\|File' | wc -l` 
  if [[ $e = 1 ]]
  then
    qsub -v NAME=$2,GWASNAME=$pheno,IDP=${idpname},IDPLIST=$4 -N ${idpname}_${pheno} run.qsub
  fi
else
  qsub -v NAME=$2,GWASNAME=$pheno,IDP=${idpname},IDPLIST=$4 -N ${3}_${pheno} run.qsub
fi
for i in $(seq 1 "${nrepeat}"); do
  idpname="${3}.n${nrepeat}_${i}"
  ff=logs/run_${2}_${pheno}_${idpname}.log
  if [[ -f $ff ]]
  then
    e=`cat $ff | tail -n 1 | grep 'failed\|kill\|Errno\|File' | wc -l` 
    if [[ $e = 1 ]]
    then
      qsub -v NAME=$2,GWASNAME=$pheno,IDP=${idpname},IDPLIST=$4,\
NREPEAT=${nrepeat},REPEAT_IDX=$i,SEED=$6,LDBLOCK=$7 \
        -N ${idpname}_${pheno} run_w_perm.qsub
    fi
  else
    qsub -v NAME=$2,GWASNAME=$pheno,IDP=${idpname},IDPLIST=$4,\
NREPEAT=${nrepeat},REPEAT_IDX=$i,SEED=$6,LDBLOCK=$7 \
      -N ${idpname}_${pheno} run_w_perm.qsub
  fi
done

