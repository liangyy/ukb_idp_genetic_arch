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
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas"
idp_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas"

for rand in "${rand_seeds[@]}"; do
  pheno_list="logs/group2.rand_${rand}.pheno_list.txt"
  ls "${gwas_dir}/group2.rand_${rand}"* | \
    sed "s#${gwas_dir}/##g" | \
    sed 's/group2.//g' | \
    sed 's/.txt.gz//g' \
    > "${pheno_list}"
  for pheno in $(cat "${pheno_list}"); do
    idp_tag="group1.rand_${rand}"
    idp_list="logs/${idp_tag}.idp_list.txt"
    if [[ -f "${idp_list}" ]]; then
      rm "${idp_list}"
    fi
    tmp=($(ls "${idp_dir}/trans_qtl.param1.group_group1.rand_${rand}.omed.chr22/"))
    for tt in "${tmp[@]}"; do
      echo "${idp_dir}/trans_qtl.param1.group_group1.rand_${rand}.omed.chr{chr_num}/${tt}" >> "${idp_list}"
    done
    fout="logs/run_${midname}_${pheno}_${idp_tag}.out"
    if [[ -f "${fout}" ]]; then
      msg="$(cat "${fout}" | grep Exit | tail -n 1 | grep 1)"
      if [[ ! -z "${msg}" ]]; then
        qsub -v \
          NAME="${midname}",\
GWASNAME="${pheno}",\
IDP="${idp_tag}",\
IDPLIST="${idp_list}",\
OUTDIR="${outdir}" \
          -N "${idp_tag}_${pheno}" \
          run.qsub
      fi
    else
      qsub -v \
        NAME="${midname}",\
GWASNAME="${pheno}",\
IDP="${idp_tag}",\
IDPLIST="${idp_list}",\
OUTDIR="${outdir}" \
        -N "${idp_tag}_${pheno}" \
        run.qsub
    fi
  done
done
