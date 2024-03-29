# args1: nametag
nametag=$1
phenopre=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/gw_lasso/dmri_4th_round/list_misc/$nametag/pheno_list
batchsize=4
filein=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/fourth_round_idp_preprocessing/$nametag.parquet

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
  if [[ -f logs_en/$nametag/run_$i.out ]]
  then
    tmp=`cat logs_en/$nametag/run_$i.log | tail -n 1 | grep 'failed\|kill\|Errno' | wc -l`
    if [[ $tmp = 1 ]]
    then
      echo qsub -v BATCH=$i,NAMETAG=$nametag -N gw_en_$i run_en.qsub
    fi
  else
    # :
    echo qsub -v BATCH=$i,NAMETAG=$nametag -N gw_en_$i run_en.qsub
  fi
done
