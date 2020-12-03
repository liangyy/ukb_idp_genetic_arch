import sys
import pandas as pd

file_in = sys.argv[1]
file_out = sys.argv[2]

df = pd.read_parquet(file_in)
df = df[ df.chr.isin(['21', '22']) ].reset_index(drop=True)
df1 = df.iloc[ :1400, : ]
df11 = df.iloc[ 2000:3000, : ]
df12 = df.iloc[ 8000:9000, : ]
df2 = df.iloc[ -3000:, : ]
df21 = df.iloc[ -5000:-3500, : ]
df22 = df.iloc[ -7000:-6500, : ]
df = pd.concat([df1, df11, df12, df2, df21, df22], axis=0).reset_index(drop=True)
df.to_parquet(file_out)

