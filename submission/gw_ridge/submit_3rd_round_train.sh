datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/third_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/gw_ridge_3rd

mkdir -p $outdir

nametag=third_round_dmri

screen -dmS $nametag bash train_generic.screen $datadir/$nametag.parquet $outdir/$nametag.gw_ridge_beta.parquet train_$nametag
