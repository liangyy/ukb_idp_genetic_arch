genotypes=(
  ["banded_200"]="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan/group1.geno_cov.chr{chr_num}.banded.npz"
  ["banded_1000"]="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan_1000/group1.geno_cov.chr{chr_num}.banded.npz"
  ["naive"]="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan/group1.geno_cov.chr{chr_num}.naive.h5"
)

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/calc_inflation"
mkdir -p "${outdir}"

for k in "${!genotypes[@]}"; do
  qsub -N "${k}_cor" \
    -v NAME="${k}_cor",\
GENO_PATTERN="${genotypes[${k}]}",\
OUTDIR="${outdir}",\
EXTRA_ARG="--correlation" \
    run.qsub
  qsub -N "${k}_cov" \
    -v NAME="${k}_cov",\
GENO_PATTERN="${genotypes[${k}]}",\
OUTDIR="${outdir}" \
    run.qsub
done
