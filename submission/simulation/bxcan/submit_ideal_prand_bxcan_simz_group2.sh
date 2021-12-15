declare -A rand_seeds
rand_seeds=(
  [1]=2
  [2]=3
  [3]=4
  [4]=5
  [5]=1)
  
nametag="ideal_prand_simz"
geno_cov="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan/group1.geno_cov.chr{chr_num}.banded.npz"

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/ideal_prand_bxcan_simz"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"

rand0=2000
kk2=0

for rand in "${!rand_seeds[@]}"; do
  rand_idp="${rand_seeds[${rand}]}"
  (( kk2 = rand0 + rand + 20 ))
  idpname="param1.rand_${rand_idp}.snp_effect"
  tagname="${nametag}_${rand}"
  for gwas in $(ls "${gwas_dir}"/group2.rand_"${kk2}".oy.* | \
    sed "s#${gwas_dir}/group2.##g" | \
    sed 's/.txt.gz//g'); do
    doit="1"
    fn="logs/run_bxcan_permz_${tagname}_${gwas}.log"
    if [[ -f "${fn}" ]]; then
      doit=""
      tmp=$(cat "${fn}" | tail -n 1 | grep 'Error\|Quit')
      if [[ ! -z "${tmp}" ]]; then
        doit="1"
      fi
    fi
    if [[ ! -z "${doit}" ]]; then
      qsub \
        -v TAGNAME="${tagname}",\
GWASNAME="${gwas}",\
GENO_COV="${geno_cov}",\
IDP_WEIGHT="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes/${idpname}.parquet",\
OUTDIR="${outdir}" \
        -N "${nametag}" \
        run_ideal_bxcan_simz.qsub
    fi
  done
done


