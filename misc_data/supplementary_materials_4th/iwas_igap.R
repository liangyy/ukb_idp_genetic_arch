# require(devtools)
# install_github("kathalexknuts/MVIWAS")
# require(MVIWAS)

outdir = 'iwas_igap'
dir.create(outdir)
# system(paste0("wget 'https://www.dropbox.com/s/m7yhcx6h6pqmlgl/0019.txt.gz?dl=0' -O ", outdir, "/0019.txt.gz"))
# system(paste0("gunzip ", outdir, "/0019.txt.gz"))
# system(paste0("wget 'https://www.dropbox.com/s/6xcofhwbnyre0s5/positions.txt.gz?dl=0' -O ", outdir, "/positions.txt.gz"))
# system(paste0("gunzip ", outdir, "/positions.txt.gz"))
# system(paste0("mv ", outdir, "/positions.txt ", outdir, '/new.positions.txt'))

system(paste0("cat ", outdir, "/0019.txt | awk -F'\t' 'NR==1{print;next}; {$4=10**(-1*$4); print}' > ", outdir, "/0019tmp.txt"))
system(paste0("paste ", outdir, "/new.positions.txt ", outdir, "/0019tmp.txt > ", outdir, "/IDP0019.txt"))
system(paste0("sed -e '1s/RSID/SNP/' -e '1s/PVAL/P/' -e '1s/CHROM/CHR/' ", outdir, "/IDP0019.txt > ", outdir, "/IDP0019tmp.txt"))

# setwd(outdir)
system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

system("plink --vcf ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --out chr21")

system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")

populations <- read.table("./integrated_call_samples_v3.20130502.ALL.panel", header = T)
EUR <- as.character(populations[populations$super_pop == "EUR","sample"])
write.table(cbind(EUR, EUR), "./EUR.txt", col.names = F, row.names = F, quote = F)

system("plink --bfile ./chr21 --keep ./EUR.txt --make-bed --out 1000G.EUR.21")

#remove duplicated variants
system("plink --bfile 1000G.EUR.21 --write-snplist --out ./all_snps")

system("cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist")

system("plink --bfile 1000G.EUR.21 --exclude duplicated_snps.snplist --make-bed --out 1000G.EUR.21.DuplicatesRemoved")