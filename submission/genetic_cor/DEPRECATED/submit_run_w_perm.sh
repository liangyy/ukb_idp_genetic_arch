# ARGS1: pheno id
# ARGS2: gwas name tag (config middle name)
# ARGS3: T1 IDP list 
# ARGS4: dMRI IDP list
# ARGS5: output dir
# ARGS6: nrepeat
# ARGS7: ldblock file path

ff="logs/run_w_perm_${2}_${1}.n${6}.err"
if [[ -f $ff ]]
then
  e=`cat $ff | tail -n 1 | grep 'failed\|kill\|Errno\|Error' | wc -l` 
  if [[ $e = 1 ]]
  then
    qsub -v \
      NAME=$2,\
GWASNAME=$1,\
IDPLIST_T1=$3,\
IDPLIST_DMRI=$4,\
NREPEAT=$6,\
OUTDIR=$5,\
LDBLOCK=$7 \
      -N $1 \
      run_w_perm.qsub
  fi
else
  qsub -v \
    NAME=$2,\
GWASNAME=$1,\
IDPLIST_T1=$3,\
IDPLIST_DMRI=$4,\
NREPEAT=$6,\
OUTDIR=$5,\
LDBLOCK=$7 \
    -N $1 \
    run_w_perm.qsub
fi

