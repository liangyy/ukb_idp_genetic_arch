mkdir -p logs
for i in `seq 1 22`
do 
  sbatch --export=CHR=$i prep_1000g_genotype.sbatch
done