outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan"
mkdir -p "${outdir}"
mkdir -p logs

for i in $(seq 1 22); do
  echo qsub -v \
    CHR="${i}",\
NAME="group1_banded",\
OUTDIR="${outdir}" \
    run_geno_cov.qsub
done