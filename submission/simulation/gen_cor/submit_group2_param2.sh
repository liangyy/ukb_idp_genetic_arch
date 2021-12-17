rand_seeds=(
  1
  2
  3
  4
  5
)
h2s=(
  0.3
  0.5
  0.7
  0.9
)

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/gen_cor_param2"
mkdir -p "${outdir}"
mkdir -p logs

midname="group2_gen_cor_param2"
gwas_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/format_oy_gwas_param2"
idp_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas_param2"

rand0=2000
kk1=0
kk2=0

for rand in "${rand_seeds[@]}"; do
  (( kk2 = rand0 + rand + 20 ))
  (( kk1 = rand0 + rand + 10 ))
  pheno_list="logs/param2.group2.rand_${kk2}.pheno_list.txt"
  ls "${gwas_dir}/group2.rand_${kk2}"* | \
    sed "s#${gwas_dir}/##g" | \
    sed 's/group2.//g' | \
    sed 's/.txt.gz//g' \
    > "${pheno_list}"
  for pheno in $(cat "${pheno_list}"); do
    h2=""
    for h2i in "${h2s[@]}"; do
      tmp=$(echo "${pheno}" | grep "h2_${h2i}")
      if [[ ! -z "${tmp}" ]]; then
        h2="${h2i}"
      fi
    done
    if [[ -z "${h2}" ]]; then
      echo "WARNING: SOMETHING WRONG FOR ${gwas}"
    fi
    idp_tag="group1.rand_${kk1}"
    idp_list="logs/param2.${idp_tag}.h2_${h2}.idp_list.txt"
    if [[ -f "${idp_list}" ]]; then
      rm "${idp_list}"
    fi
    tmp=($(ls "${idp_dir}/trans_qtl.param2.group_group1.rand_${kk1}.omed.chr22/" | grep "h2_${h2}"))
    for tt in "${tmp[@]}"; do
      echo "${idp_dir}/trans_qtl.param2.group_group1.rand_${kk1}.omed.chr{chr_num}/${tt}" >> "${idp_list}"
    done
    fout="logs/run_${midname}_${pheno}_${idp_tag}.log"
    doit="1"
    if [[ -f "${fout}" ]]; then
      doit=""
      msg="$(cat "${fout}" | tail -n 1 | grep 'Error\|lock')"
      if [[ ! -z "${msg}" ]]; then
        doit="1"
      fi
    fi
    if [[ ! -z "${doit}" ]]; then
        qsub -v \
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
