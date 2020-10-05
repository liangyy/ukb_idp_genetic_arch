phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/idp_list.txt
batchsize=2

mkdir -p logs_2
mkdir -p batch_list_2
cd batch_list_2

split -d -a 3 -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list_2/`
do
  if [[ -f logs_2/run_$i.out ]]
  then
    tmp=`cat logs_2/run_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      qsub -v BATCH=$i -N gw_lasso_$i run.qsub
    fi
  else
    :
    # qsub -v BATCH=$i -N gw_lasso_$i run.qsub
  fi
done
