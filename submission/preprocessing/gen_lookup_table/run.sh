# environment setup
source /vol/bmd/yanyul/miniconda3/etc/profile.d/conda.sh
conda activate imlabtools

infile=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/hapmap3_eur_snplist/hapmap3_eur_pop_CEU_maf_0.01_with_qc.txt
outfile=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/hapmap3_eur_snplist/hapmap3_eur_pop_CEU_maf_0.01_with_qc.lookup_table.tsv.gz 

export PYTHONPATH=/vol/bmd/yanyul/GitHub/misc-tools/liftover_snp:/vol/bmd/yanyul/GitHub/misc-tools/pyutil:$PYTHONPATH
python /vol/bmd/yanyul/GitHub/misc-tools/hapmap3_snps/gen_lookup_table.py \
  --input $infile \
  --input_map /vol/bmd/yanyul/data/hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.map \
  --liftover_chain /vol/bmd/yanyul/data/chain_files/hg18ToHg19.over.chain.gz \
  --output $outfile
