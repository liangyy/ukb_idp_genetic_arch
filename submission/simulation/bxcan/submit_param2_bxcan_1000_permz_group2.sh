rand_seeds=(
  1
  2
  3
  4
  5)
h2s=(
  0.3
  0.5
  0.7
  0.9
)
nametag="permz_1000_param2"
geno_cov="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan_1000/group1.geno_cov.chr{chr_num}.banded.npz"

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan_1000_permz_param2"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas_param2"

rand0=2000
kk1=0
kk2=0

for rand in "${rand_seeds[@]}"; do
  (( kk1 = rand0 + rand + 10 ))
  (( kk2 = rand0 + rand + 20 ))
  idpname="param2.group_group1.rand_${kk1}.ridge"
  tagname="${nametag}_${rand}"
  for gwas in $(ls "${gwas_dir}"/group2.rand_"${kk2}".oy.* | \
    sed "s#${gwas_dir}/group2.##g" | \
    sed 's/.txt.gz//g'); do
    h2=""
    for h2i in "${h2s[@]}"; do
      tmp=$(echo "${gwas}" | grep "h2_${h2i}")
      if [[ ! -z "${tmp}" ]]; then
        h2="${h2i}"
      fi
    done
    if [[ -z "${h2}" ]]; then
      echo "WARNING: SOMETHING WRONG FOR ${gwas}"
    fi
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
IDP_WEIGHT="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/train_ridge_param2/split.${idpname}.h2_${h2}.parquet",\
OUTDIR="${outdir}" \
        -N "${nametag}" \
        run_bxcan_permz.qsub
    fi
  done
done


