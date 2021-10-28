rand_seeds=(
  1
  2
  3
  4
  5)
sample_size=14409

gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas"
geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"
mkdir -p "${outdir}"

for rand in "${rand_seeds[@]}"; do
  ls "${gwas_dir}/trans_qtl.param1.group_group2.rand_${rand}.oy.chr22" | \
    sed 's/.parquet//g' > "logs/rand_${rand}.txt"
  echo qsub -v \
    PHENO_LIST="logs/rand_${rand}.txt",\
INPUT_PATTERN="${gwas_dir}/trans_qtl.param1.group_group2.rand_${rand}.oy.chr{chr_num}/{pheno}.parquet",\
SNP_BIM_PATTERN="${geno_dir}/group2.chr{chr_num}.bim",\
OUTPUT_PATTERN="${outdir}/{pheno}.txt.gz",\
SAMPLE_SIZE="${sample_size}" \
    run.qsub
done