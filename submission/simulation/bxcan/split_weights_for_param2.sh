# conda activate brainxcan

rand_seeds=(
  1
  2
  3
  4
  5)

datadir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/train_ridge_param2"  

rand0=2000
kk1=0

for rand in "${rand_seeds[@]}"; do
  (( kk1 = rand0 + rand + 10 ))
  idpname="param2.group_group1.rand_${kk1}.ridge"
  python split_weights_by_h2.py \
    --weight_prefix "${datadir}/${idpnames}" \
    --h2s 0.3 0.5 0.7 0.9 \
    --output_prefix "${datadir}/split.${idpnames}"
done
