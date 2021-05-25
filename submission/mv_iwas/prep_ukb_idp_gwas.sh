# args1: download url
# args2: filename
# args3: tag
# args4: outdir
# args5: new position file
# args6: download only

url=$1
fn=$2
tag=$3
outdir=$4
newposition=$5

cd $outdir


if [[ ! -f $fn ]]
then
  wget $url -O $fn
fi 

if [[ -z $6 && ! -f $tag.tmp2.txt ]]
then
  # gunzip $fn
  zcat $fn | awk -F'\t' 'NR==1{print;next}; {$4=10**(-1*$4); print}' > $tag.tmp.txt
  paste $newposition $tag.tmp.txt > $tag.txt
  sed -e '1s/RSID/SNP/' -e '1s/PVAL/P/' -e '1s/CHROM/CHR/' $tag.txt > $tag.tmp2.txt
fi
