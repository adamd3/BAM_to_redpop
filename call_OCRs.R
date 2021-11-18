#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 2. fseq_dir = path to directory containing fseq peak calls
## 3. genome = human genome version (either `hg19` or `hg38`)

if (!length(args)==3) {
    stop("Must provide exactly 3 arguments")
} else {
    libs_file = args[1]
    fseq_dir = args[2]
    genome = args[3]
}

if (genome=="hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome_obj = BSgenome.Hsapiens.UCSC.hg38
} else {
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome_obj = BSgenome.Hsapiens.UCSC.hg19
}

libNames <- read.table(libs_file, header=TRUE)[,1]
fseq_paths <- paste0(fseq_dir, "/", libNames, ".fseqPeaks")

libNames <- libNames[file.exists(fseq_paths)]
fseq_paths <- fseq_paths[file.exists(fseq_paths)]

call_OCRs <- function(libNames, fseq_paths, genome_obj) {
    fseqPeaks <- sapply(fseq_paths, function(x) {
        toGRanges(
            x, format = "narrowPeak",
            header = FALSE, use.names = FALSE
        )
    })
    largePeaks <- sapply(fseqPeaks, function(x) {
        x[which(width(x) > 3200)]
    })
    smallPeaks <- sapply(fseqPeaks, function(x) {
        x[which(width(x) <= 3200)]
    })
    ## set minimum peak size to 3.2 kb
    smallPeaksResized <- sapply(smallPeaks, function(x) {
        GenomicRanges::resize(
            x, width = 3200, fix = "center",
            ignore.strand = TRUE, use.names = TRUE
        )
    })
    ## combine and sort
    peaksMerged <- sapply(seq_along(smallPeaksResized), function(x) {
        sort(c(smallPeaksResized[[x]], largePeaks[[x]]))
    })
    ## merge overlapping peaks
    peaksMergedReduced <- sapply(peaksMerged, function(x) {
        GenomicRanges::reduce(x)
    })
    ## truncate coordinates at the end of each chromosome
    seqlevelsStyle(genome_obj) <- "Ensembl"
    nContigs <- length(seqnames(genome_obj))
    gGr <- GRanges(
        seqnames = seqnames(genome_obj),
        ranges = IRanges(start = rep(1, nContigs),
        end = seqlengths(genome_obj)),
        strand = rep("*", nContigs)
    )
    keepChr <- standardChromosomes(gGr)
    keepChr <- keepChr[!keepChr %in% c("chrM", "MT", "chrY", "Y")]
    gGr <- keepSeqlevels(gGr, keepChr, pruning.mode = "coarse")
    fseqPeaksProcessed <- sapply(peaksMergedReduced, function(x){
        subsetByOverlaps(x, gGr, ignore.strand = TRUE)
    })
    names(fseqPeaksProcessed) <- libNames
    fseqPeaksProcessed
}

fseqOCRs <- call_OCRs(libNames, fseq_paths, genome_obj)

saveRDS(fseqOCRs,  "fseqOCRs.rds")
