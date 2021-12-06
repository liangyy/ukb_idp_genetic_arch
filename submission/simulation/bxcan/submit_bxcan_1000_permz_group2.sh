rand_seeds=(
  1
  2
  3
  4
  5)
  
nametag="group2_bxcan_1000_permz"

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/bxcan_1000_permz"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"

for rand in "${rand_seeds[@]}"; do
  for gwas in $(ls "${gwas_dir}"/group2.rand_"${rand}".oy.* | \
    sed "s#${gwas_dir}/group2.##g" | \
    sed 's/.txt.gz//g'); do
    echo qsub \
      -v CONFIGNAME="${nametag}",\
GWASNAME="${gwas}",\
OUTDIR="${outdir}" \
      -N "group2_bxcan_${nametag}" \
      run_bxcan_permz.qsub
  done
done


