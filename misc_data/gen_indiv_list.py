'''
This script is run on Washington.
Extract the list of individual ID from the final phenotype parquet file.
The file is at 
/vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet 
which is a 24409 x 889 table with the individual IDs in column called
individual.
'''

final_pheno = '/vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet'
output = '/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/final_phenotype_indivlist/2020-05-18_final-phenotypes.indiv_list.txt'

# conda activate py3
import pandas as pd
e = pd.read_parquet(final_pheno)
indivlist = e.individual
with open(output, 'w') as f:
    for i in indivlist:
        f.write(i + '\n')

