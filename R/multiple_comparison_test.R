# The raw count data were collected using the STAR-RSEM pipeline from fastq files on GSE181961 and GSE293756.

### multiple comparison test ###
in_f  = "raw_counts_4G.txt"
out_f1 = "DEG_TCC_cal_4G.txt"

# Read the input table
master = read.table(in_f, header = TRUE, row.names = 1, sep = "\t", quote = "")

col_end   = ncol(master)
param_G1  = 3
param_G2  = 3
param_G3  = 3
param_G4  = 3
param_FDR = 0.05

ex_master = master[, 3:col_end]

data = ex_master[rowSums(ex_master) > 10, ]

# Batch-effect adjustment
library(sva)
adj_run  = c(rep(1, param_G1 + param_G2), rep(2, param_G3 + param_G4))
adjusted = ComBat_seq(as.matrix(data), batch = adj_run, group = NULL)

# Normalization and DEG with TCC
library(TCC)
data.cl = c(rep(1, param_G1),
            rep(2, param_G2),
            rep(3, param_G3),
            rep(4, param_G4))
tcc = new("TCC", adjusted, data.cl)

# Estimate normalization factors (TMM) using edgeR-based iterative procedure
tcc = calcNormFactors(tcc,
                      norm.method = "tmm",
                      test.method = "edger",
                      iteration = 3,
                      FDR = 0.1,
                      floorPDEG = 0.05)
normalized = getNormalizedData(tcc)

# Differential expression (multiple-comparisons across groups via edgeR in TCC)
tcc    = estimateDE(tcc, test.method = "edger", FDR = param_FDR)
result = getResult(tcc, sort = FALSE)

tmp = cbind(rownames(tcc$count), tcc$count, result)
ave_count = log2(rowMeans(tcc$count) + 1)
tmp$a.value = ave_count
tmp = tmp[order(tmp$rank), ]

# Write output
write.table(tmp, out_f1, sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE)