# args1: IDP file, e.g. /scratch/t.cri.yliang/ukb_idp/mv_iwas/ukb_idp_gwas/25700.tmp2.txt
# args2: chr number
# args3: plink bed prefix (ld reference), e.g. /scratch/t.cri.yliang/ukb_idp/mv_iwas/1000g/1000G.EUR.$chr.DuplicatesRemoved
# args4: outdir
# args5: output tag
# args6: pval cutoff

# module load gcc/6.2.0
# module load plink/1.90
# conda activate two_sample_mr

outdir=$4
idpfile=$1
chr=$2
bedprefix=$3
outtag=$5
pval=$6

tmpfile=$outdir/$outtag.tmp
if [[ ! -f $tmpfile ]]
then
  echo "Subsetting to chromosome $chr for $idpfile"
  awk -F'\t' -v chr="$chr" 'NR==1{print;next}$1 == chr' $idpfile > $tmpfile
fi

outfile_prefix=$outdir/$outtag
if [[ ! -f $outfile_prefix.clumped ]]
then
  plink \
    --bfile $bedprefix \
    --clump $tmpfile \
    --clump-p1 1 --clump-p2 1 \
    --clump-r2 0.10 --clump-kb 1000 \
    --out $outfile_prefix
fi

weightfile=$outdir/$outtag.CPT.pval$pval.txt
if [[ ! -f $weightfile ]]
then
  echo "Post-processing $weightfile"
  Rscript post_process_clump.R \
    --clump $outfile_prefix.clumped \
    --weight $tmpfile \
    --output $weightfile \
    --pval $pval
fi

