rand_seeds=(
  1
  2
  3
  4
  5)

geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes_param2"
mkdir -p "${outdir}"
mkdir -p logs

for rand in "${rand_seeds[@]}"; do
  # group1 is used to load the SNP list
  qsub -v \
    CONFIG_MIDNAME=param2,\
OUTDIR="${outdir}",\
RAND="${rand}",\
GENO_PATTERN="${geno_dir}/group1.chr{chr_num}" \
    run_es.qsub
done
