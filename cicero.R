library(cicero)
library(viridisLite)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(Seurat)
library(Signac)
source("get_mat.R")

set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./cicero_img/", showWarnings = FALSE)

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
# Only valid regions
row.names(mat) = sub("_", "-", sub("_", ":", row.names(mat)))
mat = mat[sapply(rownames(mat), function(x) length(strsplit(x, ":")[[1]]) == 2),]
# binarize the matrix
mat@x[mat@x > 0] = 1

# format cell info
cellinfo = data.frame(cells=colnames(mat))
row.names(cellinfo) = colnames(mat)

# format peak info
peakinfo = as.data.frame(t(sapply(strsplit(row.names(mat), "_|:|-"), head)))
colnames(peakinfo) = c("chr", "bp1", "bp2")
row.names(peakinfo) = row.names(mat)

# make CDS
write("Making CDS", stderr())
write("Making fd", stderr())
fd = methods::new("AnnotatedDataFrame", data = peakinfo)
write("Making pd", stderr())
pd = methods::new("AnnotatedDataFrame", data = cellinfo)
write("New cell data set", stderr())
input_cds =  suppressWarnings(newCellDataSet(mat,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily = "binomialff"

#Ensure there are no peaks included with zero reads
write("Removing zero read peaks", stderr())
input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# Cicero CDS creation
write("Detecting genes", stderr())
input_cds = detectGenes(input_cds)
write("Estimating sizes", stderr())
input_cds = estimateSizeFactors(input_cds)
write("Dim reduction", stderr())
input_cds = reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none",
                            check_duplicates = FALSE)

tsne_coords = t(reducedDimA(input_cds))
row.names(tsne_coords) = row.names(pData(input_cds))
write("Making Cicero CDS", stderr())
cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

# Run cicero
data("human.hg19.genome")
genome = human.hg19.genome
write("Running Cicero", stderr())
conns = run_cicero(cicero_cds, genome)

# Get gene activity scores

# Get GTF
# load in your data using rtracklayer
write("Annotating with GTF", stderr())
## gene_anno = rtracklayer::readGFF("/home/gw/research/genomes/Homo_sapiens.GRCh37.82.gtf")
## data(gene_annotation_sample)
## gene_anno = gene_annotation_sample

# download and unzip
temp = tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/homo_sapiens/Homo_sapiens.GRCh37.65.gtf.gz", temp)
gene_anno = rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match plotting requirements
gene_anno$chromosome = paste0("chr", gene_anno$seqid)
gene_anno$gene = gene_anno$gene_id
gene_anno$transcript = gene_anno$transcript_id
gene_anno$symbol = gene_anno$gene_name

#### Add a column for the fData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos = subset(gene_anno, strand == "+")
pos = pos[order(pos$start),]
pos = pos[!duplicated(pos$transcript),] # remove all but the first exons per transcript
pos$end = pos$start + 1 # make a 1 base pair marker of the TSS

neg = subset(gene_anno, strand == "-")
neg = neg[order(neg$start, decreasing = TRUE),]
neg = neg[!duplicated(neg$transcript),] # remove all but the first exons per transcript
neg$start = neg$end - 1

gene_annotation_sub = rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates
# and the gene name
## gene_annotation_sub = gene_annotation_sub[,c(1:3, 8)]
gene_annotation_sub = gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] = "gene"

write("Annotating CDS by site", stderr())
input_cds = annotate_cds_by_site(input_cds, gene_annotation_sub)

#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
write("Building activity matrix", stderr())
unnorm_ga = build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga = unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes = pData(input_cds)$num_genes_expressed
names(num_genes) = row.names(pData(input_cds))

# normalize
write("Normalizing", stderr())
# cicero_gene_activities = normalize_gene_activities(unnorm_ga, num_genes)  # Reviewer suggests removing
cicero_gene_activities = unnorm_ga

# Now switch to Seurat for analysis
write("Seurat pipeline", stderr())
mat = CreateSeuratObject(counts = cicero_gene_activities, assay = "ATAC", project = "ATAC")

# Preprocess
mat = RunTFIDF(mat)
mat = FindTopFeatures(mat, min.cutoff = "q0")
mat = RunSVD(mat)

# Dimensionality reduction
mat = RunUMAP(object = mat, reduction = 'lsi', dims = 2:30)

# Cluster
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
       file = paste0(c("./cicero_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./cicero_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
