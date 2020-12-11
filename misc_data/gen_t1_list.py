import re, yaml
import pandas as pd

output = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/t1_list.txt'
output_yaml = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/t1_list.yaml'

e = pd.read_parquet('/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_T1.parquet')
mylist = list(e.columns[1:])

dict_ = {}
with open(output, 'w') as f:
    for i in mylist:
        dict_[re.sub('-', 'x', i)] = i 
        f.write(i + '\n')
        
with open(output_yaml, 'w') as f:
    yaml.dump(dict_, f)

