rand_seeds=(
  1
  2
  3
  4
  5)

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/gen_cor"
mkdir -p "${outdir}"
mkdir -p logs

midname="group2_gen_cor"
pheno_list='logs/group2.pheno_list.txt'
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"
idp_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas"
ls "${gwas_dir}" | sed 's/group2.//g' | sed 's/.txt.gz//g' > "${pheno_list}"

for rand in "${rand_seeds[@]}"; do
  for pheno in $(cat "${pheno_list}"); do
    idp_tag="group1.rand${rand}"
    idp_list="logs/${idp_tag}.idp_list.txt"
    ls "${idp_dir}/trans_qtl.param1.group_group1.rand_${rand}.omed.chr22"* | \
      sed 's/chr22/chr{chr_num}/g' > "${idp_list}"
    echo qsub -v \
      NAME="${midname}",\
GWASNAME="${pheno}",\
IDP="${idp_tag}",\
IDPLIST="${idp_list}",\
OUTDIR="${outdir}" \
      -N "${idp_tag}_${pheno}" \
      run.qsub
  done
done
