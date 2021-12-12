decalre -A rand_seeds
rand_seeds=(
  [1]=2
  [2]=3
  [3]=4
  [4]=5
  [5]=1)
  
nametag="ideal_prand_simz"
geno_cov="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan/group1.geno_cov.chr{chr_num}.banded.npz"

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/ideal_bxcan_simz"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"

for rand in "${!rand_seeds[@]}"; do
  rand_idp="${rand_seeds[${rand}]}"
  idpname="param1.rand_${rand_idp}.snp_effect"
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
      run_ideal_bxcan_simz.qsub
  done
done


