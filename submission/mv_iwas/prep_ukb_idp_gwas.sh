# args1: download url
# args2: filename
# args3: tag
# args4: outdir
# args5: new position file

url=$1
fn=$2
tag=$3
outdir=$4
newposition=$5

cd $outdir

wget $url -O $fn
gunzip 0019.txt.gz
zcat $fn | awk -F'\t' 'NR==1{print;next}; {$4=10**(-1*$4); print}' > $tag.tmp.txt
paste $newposition $tag.tmp.txt > $tag.txt
sed -e '1s/RSID/SNP/' -e '1s/PVAL/P/' -e '1s/CHROM/CHR/' $tag.txt > $tag.tmp2.txt

