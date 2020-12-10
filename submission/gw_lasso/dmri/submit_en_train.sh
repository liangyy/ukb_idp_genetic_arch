phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.txt
batchsize=2

mkdir -p logs_en
mkdir -p batch_list_en
cd batch_list_en

split -d -a 3 -l $batchsize $phenolist batch

cd ../
nowdir=`pwd`

for i in `ls batch_list_en/`
do
  if [[ -f logs_en/run_train_$i.out ]]
  then
    tmp=`cat logs_en/run_train_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      qsub -v BATCH=$i -N gw_en_$i run_en_train.qsub
    fi
  else
    # :
    qsub -v BATCH=$i -N gw_en_$i run_en_train.qsub
  fi
done
