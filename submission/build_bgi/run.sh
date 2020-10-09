bgen_prefix=/vol/bmd/data/ukbiobank/genotypes/v3/ukb_imp_chr
bgen_suffix=_v3.bgen

bgenix=/vol/bmd/yanyul/softwares/from_bitbucket/gavinband-bgen-4e33223a8dc4/build/apps/bgenix

outdir=/vol/bmd/yanyul/UKB/ukb_imp_bgi
mkdir -p $outdir
prefix=ukb_imp_chr

for chr in `seq 1 22`
do
  ln -s $bgen_prefix$chr$bgen_suffix $outdir/$prefix$chr$bgen_suffix
  $bgenix -g $outdir/$prefix$chr$bgen_suffix -index 
  unlink $outdir/$prefix$chr$bgen_suffix
done
