Here we rely on the `ldsc.py` script in [this branch](https://github.com/liangyy/ldsc/tree/massive_cor) where we made minor changes to run with the tensorqtl parquet directly.

**Computing environment**: 

To enable loading parquet, we need `pandas.read_parquet` and it requires `pandas>=0.21`. So, we require `pandas=0.21` and accordingly `fastparquet` and `python-snappy`.
(*Book-keeping*: On CRI, use `conda activate ldsc`).
