source("get_mat.R")

args = commandArgs(TRUE)
seed = args[1]
n = args[2]
input = args[3]
output = args[4]

set.seed(seed)

# Load matrix.
mat = read10Xs(input)

# Subset matrix.
matSubset = mat[,sample(1:ncol(mat), n, replace = FALSE)]

# Save matrix.
dir.create(output, recursive = TRUE)

# expression matrix
mH = file.path(paste0(output, "/matrix.mtx"))
Matrix::writeMM(matSubset, file = mH)
R.utils::gzip(mH, overwrite = TRUE)

# data frame of genes
write.table(row.names(matSubset), file = gzfile(paste0(output, "/features.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data frame of cell barcodes
write.table(gsub("\\.", "-", colnames(matSubset)), file = gzfile(paste0(output, "/barcodes.tsv.gz")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
