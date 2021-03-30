datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/third_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability_3rd_round

mkdir -p $outdir

nametag=third_round_dmri

bash run_generic.sh $datadir/$nametag.parquet $outdir $nametag
