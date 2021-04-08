# args1: nametag
nametag=$1
outdir=/scratch/t.cri.yliang/ukb_idp/gw_elastic_net_t1_4th/$nametag
phenopre=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/gw_lasso/t1_4th_round/list_misc/$nametag/pheno_list
tmplist=batch_list_$nametag.txt
outdir2=/scratch/t.cri.yliang/ukb_idp/gw_elastic_net_t1_4th
if [[ -f $tmplist ]]
then
  rm $tmplist
fi

for l in `ls batch_list_en/$nametag/`
do
  echo "$outdir/$l.gw_lasso.weights.tsv.gz" >> $tmplist
done

# conda activate ukb_idp
export PYTHONPATH=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/pyutil

rename_yaml=$phenopre.yaml

python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/gw_lasso/post_process_training.py \
  --model_list $tmplist \
  --rename_yaml $rename_yaml \
  --output $outdir2/$nametag.gw_elastic_net_beta.parquet
  
