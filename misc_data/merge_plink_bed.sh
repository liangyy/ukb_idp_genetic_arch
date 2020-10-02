plink_exec=/vol/bmd/yanyul/softwares/plink  # plink1.9 for now since plink2 does not support --merge-list yet
plink2_exec=/vol/bmd/yanyul/softwares/plink2
memory_in_mb=60000
genotype_prefix=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.chr

# this list of snps were obtained in initial plink1.9 call which returned an error and suggested to exclude them.
multi_allelic_snplist=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all-merge.missnp

outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes
outprefix=IDP_HM3_finalPheno.merged_all
tempdir=$outdir/temp

# go to working directory
cd /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/misc_data

# make temporary folder
mkdir -p $tempdir

# step1: merge 22 chromosomes using plink1.9

tmpfile_bedlist=$tempdir/tmpfile_bedlist.txt
if [[ -f $tmpfile_bedlist ]]
then
  rm $tmpfile_bedlist
fi

for i in `seq 1 22`
do
  echo $tempdir/chr$i >> $tmpfile_bedlist 
  # need to exclude multi allelic snps
  $plink_exec \
    --bfile $genotype_prefix$i \
    --exclude $multi_allelic_snplist \
    --make-bed \
    --out $tempdir/chr$i
done

$plink_exec \
  --merge-list $tmpfile_bedlist \
  --make-bed \
  --memory $memory_in_mb \
  --out $outdir/$outprefix > merge_plink_bed_step1.log 2>&1 

# remove intermediate files
rm -r $tempdir

# step2: convert BED to PGEN
$plink2_exec \
  --bfile $outdir/$outprefix \
  --make-pgen vzs \
  --memory $memory_in_mb \
  --out $outdir/$outprefix > merge_plink_bed_step2.log 2>&1

