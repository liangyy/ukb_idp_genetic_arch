#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -e logs/merge_t1_ridge_3rd.err
#PBS -o logs/merge_t1_ridge_3rd.out


source ~/.bash_profile
source ~/.bashrc
conda activate ukb_idp

python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/prediction/post_merge.py \
  --input_pattern /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/predicted_idp/pred_idp.t1_ridge_3rd.chr{chr_num}.parquet \
  --indiv_col indiv \
  --output /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/predicted_idp/pred_idp.t1_ridge_3rd.parquet \
  > ${PBS_O_WORKDIR}/logs/merge_t1_ridge_3rd.log 2>&1
  
