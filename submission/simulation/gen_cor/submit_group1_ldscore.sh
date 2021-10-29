outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/gen_cor"
mkdir -p "${outdir}"

geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
geno_prefix="group1.chr"
geno_suffix=""

for i in $(seq 1 22); do
  qsub -v \
    CHR="${i}",\
GENO_PREFIX="${geno_dir}/${geno_prefix}",\
GENO_SUFFIX="${geno_suffix}",\
OUTDIR="${outdir}",\
NAME="group2" \
    calc_ldscore.qsub
done