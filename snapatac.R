library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(RColorBrewer)
source("get_mat.R")

set.seed(0) # Set seed for entire stochastic process.

## sink(file="/dev/null") # No other messages
dir.create("./snapatac_img/", showWarnings = FALSE)

args = commandArgs(TRUE)
labelPath = args[1]
dashFlag = args[2]
blackListFile = args[3]
output = args[4]
inputs = args[5:length(args)]

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

if (file.exists(paste0(inputs[1], "/matrix.mtx.gz"))) {
  mat = as(read10Xs(inputs), "dgCMatrix")
  peakinfo = as.data.frame(t(sapply(strsplit(row.names(mat), "_|:|-"), head)))
  colnames(peakinfo) = c("chr", "bp1", "bp2")
  peaks = GRanges(peakinfo$chr, IRanges(start=as.numeric(as.character(peakinfo$bp1)), end=as.numeric(as.character(peakinfo$bp2))))

  # Remove empty columns or diffusion maps complains
  mat = mat[,colSums(mat) != 0]
  write(toString(dim(mat)), stderr())

  x.sp = createSnapFromPmat(t(mat), barcodes = colnames(mat), peaks = peaks)
  matType = "pmat"
} else {
# Make snap for each file
  for(i in c(1:length(inputs))) {
    invisible(system("mkdir -p ./snapatac_temp", intern=TRUE))
    bedFile = paste0(c("./snapatac_temp/temp_", as.character(i), ".bed.gz"), collapse = "")
    outFile = paste0(c("./snapatac_temp/temp_", as.character(i), ".snap"), collapse = "")
    # Remove previous files
    invisible(system(paste0(c("rm -f", bedFile, outFile), collapse = " "), intern=TRUE))
    # Make the bed
    invisible(system(paste0(c("gzip -d -c ", inputs[i], "| sort -k4,4 | sed '/^$/d' | gzip -c >", bedFile), collapse = " "), intern=TRUE))
    invisible(system(paste0(c("gzip -d -c ", inputs[i], "| sort -k4,4 | sed '/^$/d' | gzip -c >", "test.bed.gz"), collapse = " "), intern=TRUE))
    # run snap files using the bed file
    invisible(system(paste0(c("snaptools snap-pre --verbose=False --input-file=", bedFile, " --output-snap=", outFile, " --genome-name=hg19 --genome-size=/home/gw/research/genomes/hg19.chrom.sizes --min-mapq=30 --min-flen=50 --max-flen=1000 --keep-chrm=TRUE --keep-single=FALSE --keep-secondary=False --overwrite=True --max-num=20000 --min-cov=500 --verbose=True"), collapse = ""), intern=TRUE))
    invisible(system(paste0(c("snaptools snap-add-bmat --verbose=False --snap-file=", outFile, " --bin-size-list=5000"), collapse = ""), intern=TRUE))
    invisible(system(paste0(c("rm", bedFile), collapse = " "), intern=TRUE))
  }

  # Load data
  x.sp = createSnap(
      file=sapply(c(1:length(inputs)), function(x) paste0(c("./snapatac_temp/temp_", as.character(x), ".snap"), collapse = ""))
      , sample=sapply(c(1:length(inputs)), function(x) paste0(c("sample_", as.character(x)), collapse = ""))
      , num.cores=3
    )
  matType = "bmat"
}
## barcodes = read.csv(
##     "atac_v1_adult_brain_fresh_5k_singlecell.csv",
##     head=TRUE
##   )

# Get info -- need singlecell.csv for this.
## barcodes = barcodes[2:nrow(barcodes),]
## promoter_ratio = (barcodes$promoter_region_fragments+1) /
##     (barcodes$passed_filters + 1)
## UMI = log(barcodes$passed_filters+1, 10)
## data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio)
## barcodes$promoter_ratio = promoter_ratio

# Filter barcodes based on info
## barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 &
##                                 promoter_ratio <= 0.6),]
## rownames(barcodes.sel) = barcodes.sel$barcode
## x.sp = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),]
## x.sp@metaData = barcodes.sel[x.sp@barcode,]

# Get bins
if(matType == "bmat") {
  x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=1)
}

# Make binary
x.sp = makeBinary(x.sp, mat=matType)

# Remove blacklisted locations
## blackListFile = "Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz"
## system(paste0(c("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/", blackListFile), collapse = ""))
system(paste0(c("wget ", blackListFile), collapse = ""))
blackListDownloaded = basename(blackListFile)
black_list = read.table(gzfile(blackListDownloaded, "rt"))
black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  )
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
if(length(idy) > 0){x.sp = x.sp[,-idy, mat=matType]}
system(paste0(c("rm ", blackListDownloaded), collapse = ""))

# Remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM",
                                             seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat=matType]}

# Remove overlapping bins with invariant features
if(matType == "bmat") {
  bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
} else {
  bin.cov = log10(Matrix::colSums(x.sp@pmat)+1)
}
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat=matType]

if(matType == "pmat") {
  # Needed as runDiffusionMaps only uses bmat
  # Remove empty rows or diffusion maps complains
  x.sp = x.sp[rowSums(x.sp@pmat) != 0, mat=matType]

  x.sp@feature = x.sp@peak
  x.sp@bmat = x.sp@pmat
}

# Dimensionality reduction (SnapATAC uses bmat only)
x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat=matType,
    num.eigs=50
  )

# Graph-based clustering
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:20,
    k=15
  )
x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=0
  )

# Visualize
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20,
  method="umap",
  seed.use=0
  );

items = if(as.numeric(dashFlag) == 1) {
          gsub("\\.", "-", x.sp@barcode)
        } else if (as.numeric(dashFlag) == 2){
          sub("\\.", "-", sub("\\.", "#", x.sp@barcode))
        } else {
          x.sp@barcode
        }

# Get plot data frame to use
labelDf = read.csv(labelPath, stringsAsFactors = FALSE)
labelDf$originalItem = labelDf$item
labelDf$item = toupper(labelDf$item)
plotDf = data.frame(item = as.character(items)
                  , cluster = x.sp@cluster
                  , dim1 = x.sp@umap[,1]
                  , dim2 = x.sp@umap[,2]
                    )
plotDf = merge(plotDf, labelDf, by = "item")
plotDf$item = plotDf$originalItem

# Plot
plotDf$group = plotDf$cluster
clusterP = viz(FALSE, "UMAP", plotDf)
plotDf$group = plotDf$label
labelP = viz(TRUE,"UMAP", plotDf)

ggsave(clusterP,
       file = paste0(c("./snapatac_img/", output, "_clusters.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)
ggsave(labelP,
       file = paste0(c("./snapatac_img/", output, "_labels.pdf"), collapse = ""),
       useDingbats = FALSE, width=25, height = 7)

## sink() # No other messages

# Save data
outDf = plotDf[,c("item", "cluster")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
