outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/geno_cov_evd"
geno_prefix="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes/group1.chr"
geno_suffix=".bed"
mkdir -p "${outdir}"
mkdir -p logs
jobtag="group1"

for i in $(seq 1 22); do
  if (( $i < 7)); then nbatch=200; fi
  if (( $i < 17 & $i > 6 )); then nbatch=100; fi
  if (( $i < 23 & $i > 16 )); then nbatch=50; fi
  fn="logs/run_geno_cov_evd.${jobtag}.${i}.log"
  doit="1"
  if [[ -f "${fn}" ]]; then
    tmp=$(cat "${fn}" | tail -n 1 | grep Error)
    doit=""
    if [[ ! -z "${tmp}" ]]; then
      doit="1"
    fi
  fi
  if [[ ! -z "${doit}" ]]; then
    echo qsub -v \
      INPUT_GENO="${geno_prefix}${i}${geno_suffix}",\
CHR="${i}",\
NBATCH="${nbatch}",\
NAME="${jobtag}",\
OUTPUT_PREFIX="${outdir}/${jobtag}.chr${i}" \
      run_geno_cov_evd.qsub
  fi
done
