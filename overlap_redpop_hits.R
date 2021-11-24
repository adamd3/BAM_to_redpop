#!/usr/bin/env Rscript

library(GenomicRanges)
library(redpop)
library(dendextend)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 2. res_dir = directory containing redpop results
## 3. outdir = directory for output files

in_df = args[1]
res_dir = args[2]
outdir = args[3]


##------------------------------------------------------------------------------
## Read and process redpop results
##------------------------------------------------------------------------------
chip_atac_df <- read.table(in_df, header=TRUE)
fnames <- paste0(chip_atac_df[,3],"_",chip_atac_df[,4])
chromosomes <- 1:22

redpop_files <- gsub(" ","",file.path(paste0(
    apply(expand.grid(fnames, chromosomes), 1, paste0, collapse="_chr"),
    ".redpop.rds")))

redpop_files <- redpop_files[file.exists(redpop_files)]

redpop_hits <- lapply(redpop_files, readRDS)

names(redpop_hits) <- gsub("\\.redpop\\.rds","",redpop_files)

# names(redpop_hits[[1]][[1]])
# [1] "res"         "h3k27ac.pb"  "atac.pb"     "sel.range"   "covs.min"
# [6] "peakiness"   "smooth_covs"
## - `res` contains the hit (GenomicRanges object)
## - `sel.range` contains the full open chromatin region that was scanned

redpop_res <- lapply(redpop_hits, function(x){
    res <- lapply(x, function(y){ y$res })
    do.call("c", res)  ## combine Granges objects
})


## Union of all hits
unionAll <- Reduce(union, unlist(redpop_res))
length(unionAll)
# [1] 113


## get indices of overlaps of each sample with the Union of all
intersectsUnion <- lapply(redpop_res, function(x) {
    na.omit(findOverlaps(unionAll, x, select=c("all")))
    ## returns indices of overlaps in unionAll
})

sapply(intersectsUnion, length)
# CF001_0_chr1  CF002_0_chr1 CF002_14_chr1  CF003_0_chr1 CF003_14_chr1
#           38            11            40            41            13
# CF004_0_chr1 CF004_14_chr1  CF005_0_chr1  CF006_0_chr1  CF007_0_chr1
#           76            32            31            35            29
# CF007_14_chr1 CF007_30_chr1 CF009_30_chr1  CF010_0_chr1 CF010_14_chr1
#           40            28            14            32            42
# CF011_0_chr1 CF011_14_chr1  HV001_0_chr1  HV005_0_chr1  HV007_0_chr1
#           29            32             8            37            21


## make a data frame of overlaps
unionDf <- data.frame(unionAll)
unionDf$id  <- paste0(unionDf$seqnames, "_", unionDf$start) ## unique ID for each hit

overlapsDf <- data.frame(matrix(
    nrow = length(unionDf$id),
    ncol = length(names(redpop_hits))
))
rownames(overlapsDf) <- unionDf$id
colnames(overlapsDf) <- names(redpop_hits)

overlapsDf[is.na(overlapsDf)] <- 0

## populate rows with `1` if there's an overlap
for (i in seq_along(overlapsDf)){
    olidx <- as.data.frame(intersectsUnion[[i]])[,1]
    overlapsDf[olidx, i] <- 1
}

saveRDS(unionAll, file = paste0(outdir,"/unionAll.rds"))
saveRDS(overlapsDf, file = paste0(outdir,"/overlapsDf.rds"))


##------------------------------------------------------------------------------
## Hierarchically cluster samples on presence/absence matrix:
##------------------------------------------------------------------------------
hclust_redpop <- hclust(dist(t(overlapsDf)))
dend_redpop <- as.dendrogram(hclust_redpop)

## colour factor for dendrogram
colpal <- vector("character", length(names(redpop_hits)))
colpal[grep("CF", names(redpop_hits))] <- "#E69F00"
colpal[grep("HV", names(redpop_hits))] <- "#56B4E9"

colpal2 <- colpal[order.dendrogram(dend_redpop)]
labels_colors(dend_redpop) <- colpal2
dend_redpop <- dend_redpop %>% set("branches_lwd", 4)


png(
    paste0(outdir, "/hclust_redpop.png"),
    width = 900, height = 800, pointsize = 20
)
par(mar = c(12,4,1,1))
plot(dend_redpop)
colored_bars(
    colpal,
    dend_redpop,
    rowLabels = "Group"
)
dev.off()
