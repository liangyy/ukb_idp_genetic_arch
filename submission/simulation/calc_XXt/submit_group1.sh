geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/calc_xxt"

mkdir -p "${outdir}"
mkdir -p logs

echo qsub -v JOBNAME=group1,\
GENO_PATTERN="${geno_dir}/group1.chr{chr_num}",\
OUTPREFIX="${outdir}/group1" \
  run.qsub
