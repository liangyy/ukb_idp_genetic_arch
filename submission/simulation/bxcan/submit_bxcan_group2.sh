rand_seeds=(
  1
  2
  3
  4
  5)

mkdir -p configs

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"

for rand in "${rand_seeds[@]}"; do
  nametag="param1.group_group1.rand_${rand}.ridge"
  gwas_list="$(pwd)/configs/gwas_list.group2_bxcan_${nametag}.txt"
  ls "${gwas_dir}"/group2.rand_"${rand}".oy.* | \
    sed "s#${gwas_dir}/group2.##g" | \
    sed 's/.txt.gz//g' \
    > "${gwas_list}"
  cat config.group2_bxcan.yaml | \
    sed "s#PLACEHOLDER#$nametag#g" \
    > "configs/config.group2_bxcan_${nametag}.yaml"
  qsub \
    -v NAME="group2_bxcan_${nametag}",\
GWAS_LIST="${gwas_list}",\
OUTDIR="${outdir}" \
    -N "group2_bxcan_${nametag}" \
    run_bxcan.qsub
done


