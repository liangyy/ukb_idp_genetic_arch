datadir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/second_round_idp_preprocessing
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability_2nd_round

# t1
idptype=t1
mode1s="scaled regress original"
mode2s="non_idp_covar all_covar"
mode3s="w_pc no_pc"
for mode1 in $mode1s
do
  for mode2 in $mode2s
  do 
    for mode3 in $mode3s
    do
      nametag=$idptype.$mode1.$mode2.$mode3
      echo $nametag
      echo bash run_generic.sh $datadir/$nametag.parquet $outdir $nametag
    done
  done
done

# dmri
idptype=dmri
mode1s="regress original"
mode2s="non_idp_covar all_covar"
mode3s="w_pc no_pc"
for mode1 in $mode1s
do
  for mode2 in $mode2s
  do 
    for mode3 in $mode3s
    do
      nametag=$idptype.$mode1.$mode2.$mode3
      echo $nametag
      echo bash run_generic.sh $datadir/$nametag.parquet $outdir $nametag
    done
  done
done