library(cisTopic)
library(viridisLite)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(Rtsne)
library(densityClust)
library(RColorBrewer)
source("get_mat.R", chdir = TRUE)

set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./cistopic_img/", showWarnings = FALSE)

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

# Clustering
write("Clustering", stderr())
set.seed(0) # Set seed for entire stochastic process.
DR = Rtsne(t(cellassign), pca=F, check_duplicates = FALSE)
DRdist = dist(DR$Y)
dclust = densityClust(DRdist,gaussian=T)
# From monocle
rho = quantile(dclust$rho, probs = 0.95)
delta = quantile(dclust$delta, probs = 0.95)
dclust = findClusters(dclust, rho = rho, delta = delta)

# Add cluster information
densityClust = dclust$clusters
densityClust = as.data.frame(densityClust)
rownames(densityClust) = cisTopicObject@cell.names
colnames(densityClust) = 'cluster'
densityClust[,1] = as.factor(densityClust[,1])
cisTopicObject = addCellMetadata(cisTopicObject, densityClust)
items = if(as.numeric(dashFlag) == 1) {
          gsub("\\.", "-", cisTopicObject@cell.names)
        } else if(as.numeric(dashFlag) == 2) {
          sub("\\.", "-", sub("\\.", "#", cisTopicObject@cell.names))
        } else {
          cisTopicObject@cell.names
        }

# Get plot data frame to use
write("Creating data frame", stderr())
labelDf = read.csv(labelPath, stringsAsFactors = FALSE)
labelDf$originalItem = labelDf$item
labelDf$item = make.names(labelDf$item)
plotDf = data.frame(item = items
                   , cluster = cisTopicObject@cell.data$cluster
                   , dim1 = cisTopicObject@dr$cell$Umap[,1]
                   , dim2 = cisTopicObject@dr$cell$Umap[,2]
                    )
plotDf = merge(plotDf, labelDf, by = "item")
plotDf$item = plotDf$originalItem

# Plot
plotDf$group = plotDf$cluster
clusterP = viz(FALSE, "UMAP", plotDf)
plotDf$group = plotDf$label
labelP = viz(TRUE,"UMAP", plotDf)

ggsave(clusterP,
       file = paste0(c("./cistopic_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./cistopic_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
