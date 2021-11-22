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


fseq_paths <- lapply(libNames, function(x){
    paste0(fseq_dir, "/",
        paste0(x, "/",
            dir(path = paste0(fseq_dir, "/", x, "/"), pattern="\\.npf$")
        )
    )
})


## remove sex chromosomes and mitochondrial genome:
fseq_paths <- lapply(fseq_paths, function(x){
    x <- x[!grepl("Y\\.npf",x)]
    x <- x[!grepl("X\\.npf",x)]
    x <- x[!grepl("MT\\.npf",x)]
    x
})


libkeep <- unlist(lapply(fseq_paths, function(x) length(x)>1))

fseq_paths <- fseq_paths[libkeep]
libNames <- libNames[libkeep]

call_OCRs <- function(libNames, fseq_paths, genome_obj) {
    res <- lapply(fseq_paths, function(lib){
        chr_names <- gsub("\\.npf","",gsub("\\/","",gsub(
            paste(libNames,collapse="|"),"",
            gsub(paste0(fseq_dir,"/"),"",lib)
        )))
        fseqPeaks <- sapply(lib, function(x) {
            toGRanges(
                x, format = "narrowPeak",
                header = FALSE, use.names = FALSE
            )
        })
        names(fseqPeaks) <- chr_names
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
        ## combine peaks and sort
        peaksMerged <- sapply(seq_along(smallPeaksResized), function(x) {
            sort(c(smallPeaksResized[[x]], largePeaks[[x]]))
        })
        names(peaksMerged) <- names(smallPeaksResized)
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
        gGr <- subset(gGr, seqnames %in% chr_names)
        peaksMergedReduced <- lapply(seq_along(peaksMergedReduced), function(x){
            chr <- names(peaksMergedReduced)[x]
            chrGr <- subset(gGr, seqnames == chr)
            subsetByOverlaps(peaksMergedReduced[[x]], chrGr, ignore.strand = TRUE)
        })
        names(peaksMergedReduced) <- names(peaksMerged)
        peaksMergedReduced
    })
    names(res) <- libNames
    res
}

fseqOCRs <- call_OCRs(libNames, fseq_paths, genome_obj)

saveRDS(fseqOCRs,  "fseqOCRs.rds")
