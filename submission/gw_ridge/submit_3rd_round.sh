datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/third_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/gw_ridge_3rd

mkdir -p $outdir

nametag=third_round_dmri

screen -dmS $nametag bash run_generic.screen $datadir/$nametag.parquet $outdir/$nametag.perf.tsv.gz run_$nametag
    
