rand_seeds=(
  1
  2
  3
  4
  5)
  
group=(1 2)
geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/simulate_phenotypes"
mkdir -p "${outdir}"
mkdir -p logs

rand0=2000
kk=0

for group in "${group[@]}"; do
  for rand in "${rand_seeds[@]}"; do
    (( kk = rand0 + rand + group * 10 ))
    qsub -v \
      CONFIG_MIDNAME=param1,\
INDIV_LIST="${geno_dir}/group${group}.chr22.fam",\
OUTDIR="${outdir}",\
RAND="${kk}",\
ES_PREFIX="${outdir}/param1.rand_${rand}",\
GENO_PATTERN="${geno_dir}/group${group}.chr{chr_num}",\
GROUP="group${group}" \
      run_ph.qsub
  done
done


