mkdir -p logs

for i in `cat trait_list.txt`
do
  qsub -v CONFIG=$i run.qsub
done
