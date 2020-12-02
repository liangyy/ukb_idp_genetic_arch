import sys
import pandas as pd

gwas = sys.argv[1]
bim = sys.argv[2]
out = sys.argv[3]

res = []
for cc in range(21, 23):
    df_gwas = pd.read_parquet(gwas.format(chr_num=cc))
    df_bim = pd.read_csv(bim.format(chr_num=cc), sep='\t', header=None, columns=['chr', 'variant_id', 'nn', 'pos', 'reference', 'alternative'])
    df_gwas = pd.merge(df_gwas, df_bim, on='variant_id', how='inner')
    res.append(df_gwas)
res = pd.concat(res, axis=0)
res.to_parquet(out)
