datadir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation"

mkdir -p "${outdir}"
target_dir="${outdir}/split_genotypes"
mkdir -p "${target_dir}"

# sample list split
group1="group1.samples"
group2="group2.samples"
if [[ ! -f "${target_dir}/${group1}" || ! -f "${target_dir}/${group2}" ]]; then
  # all chromosomes have the same list of samples
  input="IDP_HM3_finalPheno.chr1.fam"
  Rscript split_samples.R \
    --input "${datadir}/${input}" \
    --group1 "${target_dir}/${group1}" \
    --group2 "${target_dir}/${group2}"
fi

# split BED files
input_prefix="${datadir}/IDP_HM3_finalPheno.chr"
if [[ ! -f "${target_dir}/group1.chr22.bed" || ! -f "${target_dir}/group2.chr22.bed" ]]; then
  for i in $(seq 1 22); do
    plink \
      --bfile "${input_prefix}${i}" \
      --make-bed \
      --keep "${target_dir}/${group1}" \
      --out "${target_dir}/group1.chr${i}" \
      --threads 16 \
      --memory 64000
    plink \
      --bfile "${input_prefix}${i}" \
      --make-bed \
      --keep "${target_dir}/${group2}" \
      --out "${target_dir}/group2.chr${i}" \
      --threads 16 \
      --memory 64000
  done
fi

