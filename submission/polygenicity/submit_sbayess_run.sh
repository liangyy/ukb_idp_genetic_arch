# args1: phenotype list
# args2: config middle name

pheno=$1

for i in `cat $pheno`
do
  logf=logs/run_sbayess_$IDP.log 
  if [[ -f $logf ]]
  then
    e=`cat $logf|tail -n 1|grep 'Error\|shallow'`
    if [[ ! -z $e ]]
    then
      qsub -v IDP=$i,CONFIG=$2 -N sbayess_${i}_$2 run_sbayess.qsub 
    fi
  else
    qsub -v IDP=$i,CONFIG=$2 -N sbayess_${i}_$2 run_sbayess.qsub 
  fi
done
