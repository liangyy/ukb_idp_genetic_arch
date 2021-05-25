# args1: chr number
# args2: outdir
# args3: eur_list

# use plink1.9

chr=$1
outdir=$2
eur=$3
plink_opt="--memory 16000 --threads 2"
cd $outdir

raw_vcf=ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
if [[ ! -f $raw_vcf ]]
then
  echo "Generating $raw_vcf"
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$raw_vcf -O $raw_vcf
fi

add_snpid_vcf=chr$chr.vcf.gz
if [[ ! -f $add_snpid_vcf ]]
then
  echo "Adding variant ID"
  bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $raw_vcf | bgzip > $add_snpid_vcf
fi

all_bed_pre=chr$chr
if [[ ! -f $all_bed_pre.bed ]]
then
  echo "Generating $all_bed_pre"
  plink --vcf $add_snpid_vcf --out $all_bed_pre $plink_opt
fi


eur_bed_pre=1000G.EUR.$chr
if [[ ! -f $eur_bed_pre.bed ]]
then
  echo "Generating $eur_bed_pre"
  plink --bfile $all_bed_pre --keep $eur --make-bed --out $eur_bed_pre $plink_opt
fi

#remove duplicated variants
snplist=all_snps.$eur_bed_pre.snplist
if [[ ! -f $snplist ]]
then
  echo "Generating $snplist"
  plink --bfile $eur_bed_pre --write-snplist --out all_snps.$eur_bed_pre $plink_opt
fi

out_pre=1000G.EUR.$chr.DuplicatesRemoved
if [[ ! -f $out_pre.bed ]]
then
  echo "Generating $out_pre"
  cat all_snps.$eur_bed_pre.snplist | sort | uniq -d > duplicated_snps.$eur_bed_pre.snplist
  plink --bfile $eur_bed_pre --exclude duplicated_snps.$eur_bed_pre.snplist --make-bed --out $out_pre $plink_opt
fi

echo "Done"
