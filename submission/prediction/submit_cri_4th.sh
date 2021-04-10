# args1: nametag
NAME=$1
prs_parquet=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_4th/$NAME.parquet

logdir=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/prediction/logs

for i in `seq 1 22`
do 
  if [[ -f $logdir/pred_cri_${NAME}_$i.log ]]
  then
    e=`cat $logdir/pred_cri_${NAME}_$i.log|tail -n 1 | grep 'sqlite3\|ModuleNotFoundError'` 
    if [[ ! -z $e ]]
    then 
      echo $i
    fi
  else
    echo $i
  fi
done > tmp

for i in `cat tmp`
do 
  qsub -v prs_parquet=$prs_parquet,NAME=$NAME,CHR=$i -N pred_cri_${NAME}_$i pred_cri.qsub  
done
rm tmp
