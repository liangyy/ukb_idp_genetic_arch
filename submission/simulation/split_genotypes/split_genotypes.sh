datadir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation"

mkdir -p "${outdir}"

# sample list split
if [[ ! -f "${outdir}/${group1}" || ! -f "${outdir}/${group2}" ]]; then
  group1="group1.samples"
  group2="group2.samples"
  # all chromosomes have the same list of samples
  input="IDP_HM3_finalPheno.chr1.fam"
  Rscript split_samples.R \
    -i "${datadir}/${input}" \
    -g1 "${group1}" \
    -g2 "${group2}"
fi

# split BED files
target_dir="${outdir}/split_genotypes" 
input_prefix="IDP_HM3_finalPheno.chr"
if [[ ! -d "${target_dir}" ]]; then
  mkdir "${target_dir}"
  for i in $(seq 1 22); do
    echo plink \
      --bfile "${input_prefix}${i}" \
      --keep "${group1}" \
      --out "${target_dir}/group1.chr${i}"
    echo plink \
      --bfile "${input_prefix}${i}" \
      --keep "${group2}" \
      --out "${target_dir}/group2.chr${i}"
  done
fi