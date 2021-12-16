# conda activate brainxcan
# generate a table with header 
# R2	Pearson	Spearman	phenotype

rand_seeds=(
  1
  2
  3
  4
  5)

datadir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes_param2"  
prefix="param2.rand_"
suffix=".snp_effect"
  
for rand in "${rand_seeds[@]}"; do
  python gen_perf_from_snp_effects.py \
    --parquet "${datadir}/${prefix}${rand}${suffix}.parquet" \
    --output "${datadir}/${prefix}${rand}${suffix}.perf.tsv.gz"
done
