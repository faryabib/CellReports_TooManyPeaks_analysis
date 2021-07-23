set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./signac_img/", showWarnings = FALSE)

suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(cowplot))
source("get_mat.R")

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

# object setup
mat = CreateSeuratObject(counts = mat, assay = "peaks", project = "10x_ATAC")

# preprocessing
mat = RunTFIDF(mat)
mat = FindTopFeatures(mat, min.cutoff = "q0")
mat = RunSVD(mat)

# reduction
mat = RunUMAP(object = mat, reduction = 'lsi', dims = 2:30)
mat = FindNeighbors(object = mat, reduction = 'lsi', dims = 2:30)
mat = FindClusters(object = mat, verbose = FALSE, algorithm = 3)

items = if(as.numeric(dashFlag) == 1) {
          gsub("\\.", "-", colnames(mat))
        } else if(as.numeric(dashFlag) == 2) {
          sub("\\.", "-", sub("\\.", "#", colnames(mat)))
        } else {
          colnames(mat)
        }

# Get plot data frame to use
labelDf = read.csv(labelPath)
labelDf$originalItem = labelDf$item
labelDf$item = make.names(labelDf$item)

plotDf = data.frame(item = items
                  , cluster = Idents(mat)
                  , dim1 = Embeddings(mat, reduction = "umap")[,1]
                  , dim2 = Embeddings(mat, reduction = "umap")[,2]
                    )
plotDf = merge(plotDf, labelDf, by = "item")
plotDf$item = plotDf$originalItem

# Plot
plotDf$group = plotDf$cluster
clusterP = viz(FALSE, "UMAP", plotDf)
plotDf$group = plotDf$label
labelP = viz(TRUE,"UMAP", plotDf)

ggsave(clusterP,
       file = paste0(c("./signac_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./signac_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
