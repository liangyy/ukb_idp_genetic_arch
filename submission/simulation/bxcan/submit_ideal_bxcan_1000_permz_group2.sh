rand_seeds=(
  1
  2
  3
  4
  5)
  
nametag="ideal_permz_1000"
geno_cov="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan_1000/group1.geno_cov.chr{chr_num}.banded.npz"

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/ideal_bxcan_1000_permz"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"

for rand in "${rand_seeds[@]}"; do
  idpname="param1.rand_${rand}.snp_effect"
  tagname="${nametag}_${rand}"
  for gwas in $(ls "${gwas_dir}"/group2.rand_"${rand}".oy.* | \
    sed "s#${gwas_dir}/group2.##g" | \
    sed 's/.txt.gz//g'); do
    echo qsub \
      -v TAGNAME="${tagname}",\
GWASNAME="${gwas}",\
GENO_COV="${geno_cov}",\
IDP_WEIGHT="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes/${idpname}.parquet",\
OUTDIR="${outdir}" \
      -N "${nametag}" \
      run_ideal_bxcan_permz.qsub
  done
done


