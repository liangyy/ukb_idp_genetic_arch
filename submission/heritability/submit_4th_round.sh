datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/fourth_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability_4th_round

mkdir -p $outdir

tags="fourth_round.dmri_no_pc fourth_round.dmri_w_pc fourth_round.t1_no_pc fourth_round.t1_w_pc"

for tag in $tags
do
  nametag=$tag
  bash run_generic.sh $datadir/$nametag.parquet $outdir $nametag
done