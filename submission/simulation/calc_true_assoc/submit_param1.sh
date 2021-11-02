rand_seeds=(
  1
  2
  3
  4
  5)

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/calc_true_assoc"
mkdir -p "${outdir}"
mkdir -p configs

eff_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes"
b1_list="configs/param1.b1_list.txt"
b2_list="configs/param1.b2_list.txt"
geno_pattern="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan/group1.geno_cov.chr{chr_num}.banded.npz"

if [[ -f "${b1_list}" ]]; then
  rm "${b1_list}"
fi
if [[ -f "${b2_list}" ]]; then
  rm "${b2_list}"
fi

for i in $(seq 1 30); do
  echo "B_${i}" >> "${b1_list}"
done
echo "b_y_null" >> "${b2_list}"

for rand in "${rand[@]}"; do 
  echo qsub -v \
    NAME="param1_rand_${rand}",\
GENO_COVAR_PATTERN="${geno_pattern}",\
B1="${eff_dir}/param1.rand_${rand}.snp_effect.parquet",\
B2="${eff_dir}/param1.rand_${rand}.snp_effect.parquet",\
B1_LIST="$(pwd)/${b1_list}",\
B2_LIST="$(pwd)/${b2_list}",\
OUTPUT="${outdir}/param1.rand_${rand}.csv" \
  run.qsub
done
