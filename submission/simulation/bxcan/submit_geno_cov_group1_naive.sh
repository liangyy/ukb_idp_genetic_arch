outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan"
mkdir -p "${outdir}"
mkdir -p logs

for i in $(seq 1 22); do
  if (( $i < 7)); then nbatch=200; fi
  if (( $i < 17 & $i > 6 )); then nbatch=100; fi
  if (( $i < 23 & $i > 16 )); then nbatch=50; fi
  qsub -v \
    CHR="${i}",\
NBATCH="${nbatch}",\
NAME="group1_naive_f32",\
OUTDIR="${outdir}" \
    run_geno_cov.qsub
done
