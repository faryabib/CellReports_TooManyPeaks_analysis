suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(Matrix)
  library(proxy)
  library(gplots)
  library(Rtsne)
  library(densityClust)
  library(irlba)
})
source("get_mat.R")

set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./cusanovich2018_img/", showWarnings = FALSE)

args = commandArgs(TRUE)
labelPath = args[1]
dashFlag = args[2]
output = args[3]
inputs = args[4:length(args)]

viz = function(isLabel, vizType, df) {
  p = ggplot(df, aes(x = dim1, y = dim2, color = factor(group))) +
    geom_point() +
    xlab(paste(vizType, "1")) +
    ylab(paste(vizType, "2")) +
    theme_cowplot(font_size = 10) +
    theme(aspect.ratio = 1)

  if(isLabel) {
    p = p + scale_color_discrete(guide = guide_legend(title = ""))
  } else {
    p = p + scale_color_discrete(guide = guide_legend(title = ""))
  }
  return(p)
}

# Load data
mat = read10Xs(inputs)

# Filter low counts
write("Filtering", stderr())
featureCounts = rowSums(mat)
mat = mat[featureCounts >= dim(mat)[2]*0.05,]
cellCounts = colSums(mat)
mat = mat[,cellCounts >= quantile(cellCounts,probs=0.1)]
mat = mat[rowSums(mat) > 0,]

# TFIDF normalization
write("Normalizing", stderr())
nfreqs = t(t(mat) / Matrix::colSums(mat))
idf = as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
mat = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs

# TSNE
write("Dimensionality reduction", stderr())
SVDtsne = irlba(mat, 50, 50)
d_diagtsne = matrix(0, nrow=50, ncol=50)
diag(d_diagtsne) = SVDtsne$d
SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
tsneTfidf = Rtsne(SVDtsne_vd,pca=F)

# Clustering
write("Clustering", stderr())
tsneDist = dist(tsneTfidf$Y)
dclust = densityClust(tsneDist,gaussian=T)
## dclust = findClusters(dclust, rho = 50, delta = 2.5)  # Cusanovich2018 recommendation
dclust = findClusters(dclust, rho = 2, delta = 2)  # For synthetic data from vignette

# Add cluster information
densityClust = dclust$clusters
densityClust = as.data.frame(densityClust)
rownames(densityClust) = colnames(mat)
colnames(densityClust) = 'cluster'
densityClust[,1] = as.factor(densityClust[,1])
items = if(as.numeric(dashFlag) == 1) {
          gsub("\\.", "-", colnames(mat))
        } else if(as.numeric(dashFlag) == 2) {
          sub("\\.", "-", sub("\\.", "#", colnames(mat)))
        } else {
          colnames(mat)
        }

# Get plot data frame to use
write("Creating data frame", stderr())
labelDf = read.csv(labelPath, stringsAsFactors = FALSE)
labelDf$originalItem = labelDf$item
labelDf$item = make.names(labelDf$item)
plotDf = data.frame(item = items
                   , cluster = densityClust$cluster
                   , dim1 = tsneTfidf$Y[,1]
                   , dim2 = tsneTfidf$Y[,2]
                   )
plotDf = merge(plotDf, labelDf, by = "item")
plotDf$item = plotDf$originalItem
plotDf = unique(plotDf)

# Plot
write("Plotting", stderr())
plotDf$group = plotDf$cluster
clusterP = viz(FALSE, "UMAP", plotDf)
plotDf$group = plotDf$label
labelP = viz(TRUE,"UMAP", plotDf)

ggsave(clusterP,
       file = paste0(c("./cusanovich2018_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./cusanovich2018_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
