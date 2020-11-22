# ARGS1: middle name of yaml file
# will submit 1 .. 22 for 22 chrs

nametag=$1
for i in `seq 1 22`
do
  if [[ -f logs/run_${nametag}_${i}.log ]]
  then
    e=`cat logs/run_${nametag}_${i}.log | tail -n 1 | grep FileNotFoundError`
    if [[ ! -z $e ]]
    then 
      qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run.qsub
    fi
  else
    qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run.qsub
  fi
done
