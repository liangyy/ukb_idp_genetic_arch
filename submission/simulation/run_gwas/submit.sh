rand_seeds=(
  1
  2
  3
  4
  5)

outdir="/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/simulation/run_gwas"
mkdir -p configs
mkdir -p "${outdir}"

for rand in "${rand_seeds[@]}"; do
  # omed GWAS on group 1
  nametag1="param1.group_group1.rand_${rand}.omed"
  # oy GWAS on group 2
  nametag2="param1.group_group2.rand_${rand}.oy"
  
  cat config.param1.yaml | sed "s#PLACEHOLDER#${nametag1}#g" > "configs/config.${nametag1}.yaml"
  cat config.param1.yaml | sed "s#PLACEHOLDER#${nametag2}#g" > "configs/config.${nametag2}.yaml"
  
  for i in $(seq 1 22); do
    echo qsub -v CHR=$i,NAME="${nametag1}",outdir="${outdir}" -N "${i}_${nametag1}_gwas" run.qsub
    echo qsub -v CHR=$i,NAME="${nametag2}",outdir="${outdir}" -N "${i}_${nametag2}_gwas" run.qsub
  done
done

