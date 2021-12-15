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

rand0=2000
kk1=0
kk2=0

for rand in "${rand_seeds[@]}"; do
  (( kk2 = rand0 + rand + 20 ))
  (( kk1 = rand0 + rand + 10 ))
  pheno_list="logs/group2.rand_${kk2}.pheno_list.txt"
  ls "${gwas_dir}/group2.rand_${kk2}"* | \
    sed "s#${gwas_dir}/##g" | \
    sed 's/group2.//g' | \
    sed 's/.txt.gz//g' \
    > "${pheno_list}"
  for pheno in $(cat "${pheno_list}"); do
    idp_tag="group1.rand_${kk1}"
    idp_list="logs/${idp_tag}.idp_list.txt"
    if [[ -f "${idp_list}" ]]; then
      rm "${idp_list}"
    fi
    tmp=($(ls "${idp_dir}/trans_qtl.param1.group_group1.rand_${kk1}.omed.chr22/"))
    for tt in "${tmp[@]}"; do
      echo "${idp_dir}/trans_qtl.param1.group_group1.rand_${kk1}.omed.chr{chr_num}/${tt}" >> "${idp_list}"
    done
    fout="logs/run_${midname}_${pheno}_${idp_tag}.log"
    doit="1"
    if [[ -f "${fout}" ]]; then
      doit=""
      msg="$(cat "${fout}" | tail -n 1 | grep Error)"
      if [[ ! -z "${msg}" ]]; then
        doit="1"
      fi
    fi
    if [[ ! -z "${doit}" ]]; then
        echo qsub -v \
          NAME="${midname}",\
GWASNAME="${pheno}",\
IDP="${idp_tag}",\
IDPLIST="${idp_list}",\
OUTDIR="${outdir}" \
          -N "${idp_tag}_${pheno}" \
          run.qsub
    fi
#     else
#       echo qsub -v \
#         NAME="${midname}",\
# GWASNAME="${pheno}",\
# IDP="${idp_tag}",\
# IDPLIST="${idp_list}",\
# OUTDIR="${outdir}" \
#         -N "${idp_tag}_${pheno}" \
#         run.qsub
#     fi
  done
done
