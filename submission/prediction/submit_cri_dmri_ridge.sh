prs_parquet=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/gw_ridge/gw_ridge_beta.dmri.default_theta_g_fold_5_5.parquet
NAME=dmri_ridge

for i in `seq 1 22`; do e=`cat ~/pred_cri_dmri_ridge_$i.log|tail -n 1 | grep sqlite3`; if [[ ! -z $e ]]; then echo $i; fi; done > tmp

for i in `cat tmp`
do 
  echo qsub -v prs_parquet=$prs_parquet,NAME=$NAME,CHR=$i -N $i.pred_ridge_dmri pred_cri.qsub  
done
rm tmp
