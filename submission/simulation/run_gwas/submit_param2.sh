rand_seeds=(
  1
  2
  3
  4
  5
)

geno_dir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/split_genotypes"
outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas_param2"
geno1="${geno_dir}/group1.chr{chr_num}"
geno2="${geno_dir}/group2.chr{chr_num}"

mkdir -p configs
mkdir -p "${outdir}"

rand0=2000
kk=0

for rand in "${rand_seeds[@]}"; do
  (( kk1 = rand0 + rand + 10 ))  
  (( kk2 = rand0 + rand + 20 ))
  # omed GWAS on group 1
  nametag1="param2.group_group1.rand_${kk1}.omed"
  # oy GWAS on group 2
  nametag2="param2.group_group2.rand_${kk2}.oy"
 
  cat config.param2.yaml | sed "s#PLACEHOLDER#${nametag1}#g" > "configs/config.${nametag1}.yaml"
  cat config.param2.yaml | sed "s#PLACEHOLDER#${nametag2}#g" > "configs/config.${nametag2}.yaml"
  
  for i in $(seq 1 22); do
    if [[ -f "logs/${nametag1}.${i}.log" ]]; then
      kk="$(cat "logs/${nametag1}.${i}.log" | tail -n 1 | grep 'Error\|lock')"
      if [[ ! -z "${kk}" ]]; then
        qsub -v CHR=$i,NAME="${nametag1}",OUTDIR="${outdir}",GENO_PATTERN="${geno1}.bed" -N "${i}_${nametag1}_gwas" run.qsub
      fi
    else
      qsub -v CHR=$i,NAME="${nametag1}",OUTDIR="${outdir}",GENO_PATTERN="${geno1}.bed" -N "${i}_${nametag1}_gwas" run.qsub
    fi
    if [[ -f "logs/${nametag2}.${i}.log" ]]; then
      kk="$(cat "logs/${nametag2}.${i}.log" | tail -n 1 | grep 'Error\|lock')"
      if [[ ! -z "${kk}" ]]; then
        echo qsub -v CHR=$i,NAME="${nametag2}",OUTDIR="${outdir}",GENO_PATTERN="${geno2}.bed" -N "${i}_${nametag2}_gwas" run.qsub
      fi
    else
      echo qsub -v CHR=$i,NAME="${nametag2}",OUTDIR="${outdir}",GENO_PATTERN="${geno2}.bed" -N "${i}_${nametag2}_gwas" run.qsub
    fi
  done
done

