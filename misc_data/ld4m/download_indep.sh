cd /project2/haky/yanyul/data/SLD4M/external_gwas
for i in `cat /project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/ld4m/ld4m_indep_gwas.txt`
do
  wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/independent_sumstats/$i.sumstats.gz
done

