import pandas as pd

output = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.txt'

e = pd.read_parquet('/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_dMRI.parquet')
mylist = list(e.columns[1:])

with open(output, 'w') as f:
    for i in mylist:
        f.write(i + '\n')
