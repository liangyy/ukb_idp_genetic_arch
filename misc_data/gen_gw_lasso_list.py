
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='gen_gw_lasso_list.py', description='''
        Generate the phenotype list and YAMl for gw_lasso.
    ''')
    parser.add_argument('--input', help='''
        The phenotype matrix in parquet format.
    ''')
    parser.add_argument('--output_prefix', help='''
        Prefix of output files.
    ''')
    args = parser.parse_args()

    import re, yaml
    import pandas as pd

    # output = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.txt'
    # output_yaml = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/dmri_list.yaml'
    # input = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_dMRI.parquet'
    
    e = pd.read_parquet(args.input)
    mylist = list(e.columns[1:])

    dict_ = {}
    with open(args.output_prefix + '.txt', 'w') as f:
        for i in mylist:
            dict_[re.sub('-', 'x', i)] = i 
            f.write(i + '\n')
        
    with open(args.output_prefix + '.yaml', 'w') as f:
        yaml.dump(dict_, f)
    
