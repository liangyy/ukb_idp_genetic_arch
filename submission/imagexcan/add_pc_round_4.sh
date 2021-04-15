conda activate ukb_idp
datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/prediction 
models="ridge elastic_net"
idps="t1 dmri"
for model in $models
do
  for idp in $idps
  do
    echo $model $idp
    python add_pc_pred_expr.py \
      --input $datadir/pred_idp.fourth_round.${idp}_no_pc.gw_${model}_beta.parquet \
      --input_pc $datadir/pred_idp.fourth_round.${idp}_w_pc.gw_${model}_beta.parquet \
      --output $datadir/pred_idp.fourth_round.${idp}.original.gw_${model}_beta.parquet
  done
done