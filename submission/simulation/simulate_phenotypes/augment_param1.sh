# conda activate brainxcan
rand_seeds=(
  1
  2
  3
  4
  5)
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes"
configname=param1

for rand in "${rand_seeds[@]}"; do
  python augment_snp_effects.py \
    --snp_effect_parquet "${outdir}/${configname}.rand_${rand}.snp_effect.parquet" \
    --rand_seed 2021 \
    --augment_size 300 \
    --output "${outdir}/${configname}.rand_${rand}.snp_effect.augmented.parquet"
done