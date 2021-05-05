library(reticulate)
use_condaenv('SPrediXcan2PTRS')
py = import('pyarrow')
pandas = import('pandas')
read_parquet = function(fn) {
  pandas$read_parquet(fn)
}
write_parquet = function(dd, fn) {
  py$parquet$write_table(py$Table$from_pandas(dd), fn)
}