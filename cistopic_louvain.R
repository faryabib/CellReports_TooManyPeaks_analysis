library(cisTopic)
library(viridisLite)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(Rtsne)
library(densityClust)
library(RColorBrewer)
library(Seurat)
library(Signac)
source("get_mat.R")

set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./cistopic_louvain_img/", showWarnings = FALSE)

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
    pal = colorRampPalette(brewer.pal(9, "Set1"))
    p = p + scale_color_manual(guide = guide_legend(title = ""), values = pal(nlevels(as.factor(df$group))))
  } else {
    p = p + scale_color_discrete(guide = guide_legend(title = ""))
  }
  return(p)
}

mat = read10Xs(inputs)
row.names(mat) = sub("_", "-", sub("_", ":", row.names(mat)))
# Only valid regions
mat = mat[sapply(rownames(mat), function(x) length(strsplit(x, ":")[[1]]) == 2),]
write("Finished matrix", stderr())

# Make cisTopic object
write("Creating cisTopic object", stderr())
cisTopicObject = createcisTopicObject(mat)

# Build models
write("Run models", stderr())
cisTopicObject = runWarpLDAModels(cisTopicObject, topic=c(30, 35, 40), seed=0, nCores=48, iterations = 500, addModels=FALSE)  # Set by reviewer

# Select model
write("Select model", stderr())
cisTopicObject = selectModel(cisTopicObject)

# Reduction
write("Run UMAP", stderr())
set.seed(0) # Set seed for entire stochastic process.
cisTopicObject = runUmap(cisTopicObject, target='cell', seed=0, method='Probability')

# Normalized assignments
write("Normalize assignments", stderr())
cellassign = modelMatSelection(cisTopicObject, 'cell', 'Probability')

# Now switch to Seurat for analysis
write("Seurat pipeline", stderr())
mat = CreateSeuratObject(counts = cellassign, assay = "ATAC", project = "ATAC")
mat[["cistopic"]] = CreateDimReducObject(embeddings = t(cellassign), key = "Topic", assay = DefaultAssay(mat))

# Preprocess
mat = FindTopFeatures(mat, min.cutoff = "q0")
mat = RunSVD(mat) # Depending on topics

# Dimensionality reduction
mat = RunUMAP(object = mat, reduction = 'lsi', dims = 2:30) # Depending on topics

# Cluster
mat = FindNeighbors(object = mat, reduction = 'lsi', dims = 2:30) # Depending on topics
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

# For plotting
mat = RunUMAP(object = mat, reduction = 'cistopic', dims = 2:30)

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
       file = paste0(c("./cistopic_louvain_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./cistopic_louvain_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
