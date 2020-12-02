e = read.table('test_files/output_test_training.sscore')
f = readRDS('test_files/output_test_training.y_insample.rds')
d = data.frame(y0 = f$QPHE$ypred, y1 = e$V3)
message('cor( y_pred from plink2, y_pred from snpnet ) = ', cor(d$y0, d$y1))
