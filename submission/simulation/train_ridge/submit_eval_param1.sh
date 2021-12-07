rand_seeds=(
  1
  2
  3
  4
  5
)
groups=(1 2)
param_tag="param1"
pheno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes"
geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/train_ridge"
mkdir -p "${outdir}"

for group in "${groups[@]}"; do
  for rand in "${rand_seeds[@]}"; do 
    geno_pattern="${geno_dir}/group${group}.chr{chr_num}"
    tag="${param_tag}.group_group${group}.rand_${rand}"
    pheno_parquet="${pheno_dir}/${tag}.omed.parquet"
    qsub -v \
      GENO_PATTERN="${geno_pattern}",\
PHENO_PARQUET="${pheno_parquet}",\
TAG="${tag}",\
CACHE_TAG="${geno_dir}/group${group}",\
OUTDIR="${outdir}" \
      eval.qsub
  done
done
