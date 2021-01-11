# on CRI
# require pandas=0.21, fastparquet, and python-snappy
conda activate ldsc

ldscdir=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ldsc
gwas=/gpfs/data/im-lab/nas40t2/Data/SummaryResults/imputed_gwas_hg38_1.1/imputed_pgc.scz2.txt.gz
snpmeta=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all.bim
qtldir=/scratch/t.cri.yliang/ukb_idp/idp_gwas_2nd/trans_qtl.dmri.original.all_covar.w_pc.chr{chr_num}
scoredir=/gpfs/data/im-lab/nas40t2/yanyul/data/eur_w_ld_chr
outdir=test_run_out

mkdir -p $outdir

if [[ ! -f $outdir/gwas_formatted.sumstats.gz ]]
then
  python $ldscdir/munge_sumstats.py \
    --sumstats $gwas \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --p pvalue \
    --signed-sumstats zscore,0 \
    --snp variant_id \
    --N-col sample_size \
    --out $outdir/gwas_formatted \
    --merge-alleles $scoredir/w_hm3.snplist
fi

echo python $ldscdir/ldsc.py \
  --rg $outdir/gwas_formatted.sumstats.gz,$qtldir/IDP-25361.parquet,$qtldir/IDP-25203.parquet \
  --ref-ld-chr $scoredir/ \
  --w-ld-chr $scoredir/ \
  --out $outdir/rg_out \
  --load-parquet \
  --snp-meta $snpmeta \
  --sample-size 24409
  
  
