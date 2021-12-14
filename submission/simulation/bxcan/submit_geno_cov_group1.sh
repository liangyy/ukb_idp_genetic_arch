outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan"
mkdir -p "${outdir}"
mkdir -p logs

for i in $(seq 1 22); do
  if (( $i < 7)); then nbatch=200; fi
  if (( $i < 17 & $i > 6 )); then nbatch=100; fi
  if (( $i < 23 & $i > 16 )); then nbatch=50; fi
  fn="logs/run_geno_cov.group1_banded.${i}.log"
  doit="1"
  if [[ -f "${fn}" ]]; then
    tmp=$(cat "${fn}" | tail -n 1 | grep Error)
    doit=""
    if [[ ! -z "${tmp}" ]]; then
      doit="1"
    fi 
  fi
  if [[ ! -z "${doit}" ]]; then
    qsub -v \
      CHR="${i}",\
NBATCH="${nbatch}",\
NAME="group1_banded",\
OUTDIR="${outdir}" \
      run_geno_cov.qsub
  fi
done
