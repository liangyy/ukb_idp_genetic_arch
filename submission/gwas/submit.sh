# ARGS1: middle name of yaml file
# will submit 1 .. 22 for 22 chrs

nametag=$1
for i in `seq 1 22`
do
  qsub -v CHR=$i -N ${i}_${nametag}_gwas run.qsub
done