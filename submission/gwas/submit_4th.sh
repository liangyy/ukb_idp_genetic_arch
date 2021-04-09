# ARGS1: nametag

nametag=$1

mkdir -p configs

cat config.4th.yaml | sed "s#PLACEHOLDER#$nametag#g" > configs/config.$nametag.yaml

for i in `seq 1 22`
do
  if [[ -f logs/run_2nd_${nametag}_${i}.log ]]
  then
    e=`cat logs/run_2nd_${nametag}_${i}.log | tail -n 1 | grep 'FileExistsError\|FileNotFoundError\|unlock'`
    if [[ ! -z $e ]]
    then 
      echo qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run_2nd.qsub
    fi
  else
    # :
    echo qsub -v CHR=$i,NAME=$nametag -N ${i}_${nametag}_gwas run_2nd.qsub
  fi
done
