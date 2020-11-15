prs_parquet=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/gw_ridge/gw_ridge_beta.t1.default_theta_g_fold_5_5.parquet
NAME=t1_ridge
logdir=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/prediction/logs

for i in `seq 1 22`; do e=`cat $logdir/pred_cri_t1_ridge_$i.log|tail -n 1 | grep sqlite3`; if [[ ! -z $e ]]; then echo $i; fi; done > tmp

for i in `cat tmp`
do 
  echo qsub -v prs_parquet=$prs_parquet,NAME=$NAME,CHR=$i -N $i.pred_ridge_t1 pred_cri.qsub  
done
rm tmp
