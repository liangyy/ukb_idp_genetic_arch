import pandas as pd

file_cov = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/phenotypes/imagexcan_covariate_round_1.parquet'
file_indiv = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/fourth_round_idp_preprocessing/fourth_round.dmri_no_pc.parquet'
outdir = '/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/additional_tables'

# load covariates
df_cov = pd.read_parquet(file_cov)
# load indiv
df_indiv = pd.read_parquet(file_indiv)

# sex coding: male - 1; female - 0
df_sub = df_cov[ df_cov.eid.isin(df_indiv.individual) ][['sex', 'age']]

# aggregate results
def count_zeros(x): return (x == 0).sum()
out = df_sub.agg({'age': ['mean', 'std'], 'sex': ['mean', 'sum', count_zeros]})
out.to_csv(f'{outdir}/demographic.csv', index=False)
