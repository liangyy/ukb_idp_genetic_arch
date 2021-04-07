datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/fourth_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/gw_ridge_4th

mkdir -p $outdir

tags="fourth_round.dmri_no_pc fourth_round.dmri_w_pc fourth_round.t1_no_pc fourth_round.t1_w_pc"
    
for nametag in $tags
do
  screen -dmS $nametag bash run_generic.screen $datadir/$nametag.parquet $outdir/$nametag.perf.tsv.gz run_$nametag
done  