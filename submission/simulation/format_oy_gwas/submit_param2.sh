rand_seeds=(
  1
  2
  3
  4
  5)
sample_size=14409

gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas_param2"
geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas_param2"
mkdir -p "${outdir}"
mkdir -p logs

rand0=2000
kk=0

for rand in "${rand_seeds[@]}"; do
  (( kk2 = rand0 + rand + 20 ))  
  ls "${gwas_dir}/trans_qtl.param2.group_group2.rand_${kk2}.oy.chr22" | \
    sed 's/.parquet//g' > "logs/rand_${kk2}.txt"
  echo qsub -v \
    PHENO_LIST="logs/rand_${kk2}.txt",\
NAME="rand_${kk2}",\
INPUT_PATTERN="${gwas_dir}/trans_qtl.param2.group_group2.rand_${kk2}.oy.chr{chr_num}/{pheno}.parquet",\
SNP_BIM_PATTERN="${geno_dir}/group2.chr{chr_num}.bim",\
OUTPUT_PATTERN="${outdir}/group2.rand_${kk2}.oy.{pheno}.txt.gz",\
SAMPLE_SIZE="${sample_size}" \
    run.qsub
done
