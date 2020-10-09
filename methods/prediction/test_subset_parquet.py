import argparse
parser = argparse.ArgumentParser(prog='test_subset_parquet.py', description='''
    Extract the first 50 SNPs per chromosome.
    And extract the first 3 traits.
''')
parser.add_argument('--input', help='''
    Input beta parquet.
''')
parser.add_argument('--output_prefix', help='''
    Output prefix of parquet and txt file.
''')
args = parser.parse_args()

import pandas as pd
import numpy as np
np.random.seed(10)

df = pd.read_parquet(args.input)
chrs = [1, 22]

out = []
for i in chrs:
    df_sub = df[df.chr == str(i)].reset_index(drop=True)
    df_sub = df_sub.iloc[np.unique(np.random.randint(0, 7000, size=(50))), :].reset_index(drop=True)
    df_sub = df_sub[df_sub.columns[:7]]
    out.append(df_sub)
out = pd.concat(out, axis=0)
out.to_parquet(args.output_prefix + '.parquet')
out.to_csv(args.output_prefix + '.txt', index=False, sep='\t')
