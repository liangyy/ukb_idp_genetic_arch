# args1: nametag
nametag=$1
phenopre=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/gw_lasso/dmri_3rd_round/list_misc/$nametag/pheno_list
batchsize=2
filein=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/third_round_idp_preprocessing/$nametag.parquet

mkdir -p logs_en
mkdir -p batch_list_en
mkdir -p list_misc
mkdir -p configs

mkdir -p logs_en/$nametag
mkdir -p batch_list_en/$nametag
mkdir -p list_misc/$nametag

if [[ ! -f $phenopre.txt ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/gen_gw_lasso_list.py --input $filein --output_prefix $phenopre
fi

cat config.en.yaml | sed "s#PLACEHOLDER#$nametag#g" > configs/config.$nametag.yaml

cd batch_list_en/$nametag

split -d -a 3 -l $batchsize $phenopre.txt batch

cd ../../
nowdir=`pwd`

for i in `ls batch_list_en/$nametag`
do
  if [[ -f logs_en/$nametag/run_train_$i.log ]]
  then
    tmp=`cat logs_en/$nametag/run_train_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      qsub -v BATCH=$i,NAMETAG=$nametag -N train_gw_en_$i run_en_train.qsub
    fi
  else
    # :
    qsub -v BATCH=$i,NAMETAG=$nametag -N train_gw_en_$i run_en_train.qsub
  fi
  # echo qsub -v BATCH=$i,NAMETAG=$nametag -N train_gw_en_$i run_en_train.qsub
done
