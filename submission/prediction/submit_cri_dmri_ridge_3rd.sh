prs_parquet=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/idp_models_3rd/third_round_dmri.gw_ridge_beta.parquet
NAME=dmri_ridge_3rd

# for i in `seq 1 22`; do e=`cat ~/pred_cri_dmri_ridge_3rd_$i.log|tail -n 1 | grep sqlite3`; if [[ ! -z $e ]]; then echo $i; fi; done > tmp

for i in `seq 1 22`  #`cat tmp`
do 
  echo qsub -v prs_parquet=$prs_parquet,NAME=$NAME,CHR=$i -N $i.pred_ridge_dmri pred_cri.qsub  
done
# rm tmp
