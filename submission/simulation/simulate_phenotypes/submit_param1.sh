rand_seeds=(
  1
  2
  3
  4
  5)
groups=(1 2)

geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes"
mkdir -p "${outdir}"

for group in "${groups[@]}"; do
  for rand in "${rand_seed[@]}"; do
    echo qsub -v \
      CONFIG_MIDNAME=param1 \
      INDIV_LIST="${geno_dir}/group${group}.chr22.fam" \
      OUTDIR="${outdir}" \
      RAND="${rand}" \
      GENO_PATTERN="${geno_dir}/group${group}.chr{chr_num}" \
      GROUP="group${group}" \
      run.qsub
  done
done
