import sys
import pandas as pd

file_in = sys.argv[1]
file_out = sys.argv[2]

df = pd.read_parquet(file_in)
df = df[ df.chr.isin(['21', '22']) ].reset_index(drop=True)

df.to_parquet(file_out)

