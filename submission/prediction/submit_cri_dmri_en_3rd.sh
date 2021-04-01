prs_parquet=/scratch/t.cri.yliang/ukb_idp/idp_models_3rd/third_round_dmri.gw_elastic_net_beta.parquet
NAME=dmri_en_3rd
logdir=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/prediction/logs

for i in `seq 1 22`
do 
  if [[ -f $logdir/pred_cri_dmri_en_3rd_$i.log ]]
  then
    e=`cat $logdir/pred_cri_dmri_en_3rd_$i.log|tail -n 1 | grep 'sqlite3\|ModuleNotFoundError'` 
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
  qsub -v prs_parquet=$prs_parquet,NAME=$NAME,CHR=$i -N $i.pred_en_dmri pred_cri.qsub  
done
rm tmp
