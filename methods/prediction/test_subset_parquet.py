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

BASE_PAIR = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def is_ambi(a, b):
    o = []
    for i, j, in zip(a, b):
        if i == BASE_PAIR[j]:
            o.append(True)
        else:
            o.append(False)
    return o

df = pd.read_parquet(args.input)
chrs = [ i for i in range(8, 13) ]

out = []
for i in chrs:
    df_sub = df[df.chr == str(i)].reset_index(drop=True)
    df_sub = df_sub.iloc[np.unique(np.random.randint(0, 7000, size=(135))), :].reset_index(drop=True)
    bad_ind = is_ambi(df_sub.a0.tolist(), df_sub.a1.tolist())
    # breakpoint()
    df_sub = df_sub[ ~np.array(bad_ind) ].reset_index(drop=True)
    df_sub = df_sub[df_sub.columns[:7]]
    out.append(df_sub)
out = pd.concat(out, axis=0)
out.to_parquet(args.output_prefix + '.parquet')
out.to_csv(args.output_prefix + '.txt', index=False, sep='\t')

