outdir=/scratch/t.cri.yliang/ukb_idp/gw_elastic_net_dmri
tmplist=batch_list.txt
if [[ -f $tmplist ]]
then
  rm $tmplist
fi

for l in `ls batch_list_en/`
do
  echo "$outdir/$l.gw_lasso.weights.tsv.gz" >> $tmplist
done

# conda activate ukb_idp
export PYTHONPATH=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/pyutil

rename_yaml=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.yaml

python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/gw_lasso/post_process_training.py \
  --model_list $tmplist \
  --rename_yaml $rename_yaml \
  --output $outdir/gw_lasso_beta.parquet
  
