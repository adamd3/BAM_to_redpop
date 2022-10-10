#!/usr/bin/env Rscript

library(GenomicRanges)
library(redpop)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 2. atacDatDir = Directory containing ATAC-Seq BigWig files
## 3. chipDatDir = Directory containing ChIP-Seq BigWig files
## 4. bwExt = File extension for BigWigs
## 5. chr_select = Selected chromosome to scan (as a single number) [`all` for all chromosomes]



in_df = args[1]
atacDatDir = args[2]
chipDatDir = args[3]
bwExt = args[4]
chr_select = args[5]


## minor modifications to the original functions from Ernest Turro to avoid
## conflicts between function names in different packages:
source("redpop_v2.R")
source("redpop_helpers.R")



##------------------------------------------------------------------------------
## read in data
##------------------------------------------------------------------------------
fseqOCRs <- readRDS("fseqOCRs.rds")

chip_atac_df <- read.table(in_df, header=TRUE)

### only scan available files
chip_atac_df <- subset(chip_atac_df,
    ATAC_sample %in% names(fseqOCRs) | Chip_sample %in% names(fseqOCRs)
)

atacF <- paste0(atacDatDir, "/", chip_atac_df[["ATAC_sample"]], bwExt)
chipF <- paste0(chipDatDir, "/", chip_atac_df[["Chip_sample"]], bwExt)


##------------------------------------------------------------------------------
## Run RedPop
##------------------------------------------------------------------------------
lapply(1:nrow(chip_atac_df), function(x){
    ## subset to OCRs in this sample (one per chromosome)
    sample_OCRs <- fseqOCRs[chip_atac_df[x,][["ATAC_sample"]]][[1]]
    if(!chr_select=="all"){
        sample_OCRs <- subset(sample_OCRs, names(sample_OCRs)==chr_select)
    }
    atac_file <- atacF[x]
    chip_file <- chipF[x]
    ## idenitifier for the sample
    sample <- paste0(chip_atac_df[,3],"_",chip_atac_df[,3])[x]
    if(length(sample_OCRs)>1){
        lapply(seq_along(sample_OCRs), function(chr_idx){
            chr_name <- names(sample_OCRs)[chr_idx]
            outf <- paste0(sample, "_chr", chr_name)
            chr_ocrs <- sample_OCRs[chr_idx][[1]]
            chr_res <- lapply(seq_along(chr_ocrs), function(y){
                tmp <- tryCatch(
                    redpop_v2(atac_file, chip_file, chr_ocrs[y]),
                    error = function(e) e
                )
                if(!inherits(tmp, "error")){
                    res <- NA
                    if (length(tmp$res) > 0) res <- tmp
                } else {
                    res <- NA
                }
                return(res)
            })
            chr_res <- chr_res[!is.na(chr_res)]
            saveRDS(chr_res, file = paste0(outf,".redpop.rds"))
        })
    } else {
        chr_name <- chr_select
        outf <- paste0(sample, "_chr", chr_name)
        chr_ocrs <- sample_OCRs[[1]]
        chr_res <- lapply(seq_along(chr_ocrs), function(y){
            tmp <- tryCatch(
                redpop_v2(atac_file, chip_file, chr_ocrs[y]),
                error = function(e) e
            )
            if(!inherits(tmp, "error")){
                res <- NA
                if (length(tmp$res) > 0) res <- tmp
            } else {
                res <- NA
            }
            return(res)
        })
        chr_res <- chr_res[!is.na(chr_res)]
        saveRDS(chr_res, file = paste0(outf,".redpop.rds"))
    }
})


## Scan just the AURKAIP hit in CF002:

#                ATAC_sample            Chip_sample Individual Time
# 2 SLX16960CF0020MonoATAC_A SLX20032i704i506Chip_A      CF002    0
# 3  SLX16960Ad002MonoATAC_A SLX20032i710i507Chip_A      CF002   14

atac_samps <- subset(chip_atac_df, Individual=="CF002")[["ATAC_sample"]]
chip_samps <- subset(chip_atac_df, Individual=="CF002")[["Chip_sample"]]

atacF <- paste0(atacDatDir, "/", atac_samps, bwExt)
chipF <- paste0(chipDatDir, "/",chip_samps, bwExt)

sample_OCRs_Day0 <- (fseqOCRs[atac_samps[1]][[1]])[1][[1]]
sample_OCRs_Day14 <- (fseqOCRs[atac_samps[2]][[1]])[1][[1]]

length(sample_OCRs_Day0)
# [1] 12575
length(sample_OCRs_Day14)
# [1] 15892


## get the index of the selected hit:
grtest <- GRanges(
    seqnames = "1", ranges=IRanges(start=1273785,width=1),strand = c('*'))

selected_hit <- as.data.frame(findOverlaps(sample_OCRs_Day0,grtest))[,1]
selected_hit
# [1] 41

sample_OCRs_Day0[41]
# GRanges object with 1 range and 0 metadata columns:
#       seqnames          ranges strand
#          <Rle>       <IRanges>  <Rle>
#   [1]        1 1271713-1275426      *
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths


## get the index of the selected hit:
grtest <- GRanges(
    seqnames = "1", ranges=IRanges(start=1273785,width=1),strand = c('*'))

selected_hit <- as.data.frame(findOverlaps(sample_OCRs_Day14,grtest))[,1]
selected_hit
# [1] 46

sample_OCRs_Day14[46]
# GRanges object with 1 range and 0 metadata columns:
#       seqnames          ranges strand
#          <Rle>       <IRanges>  <Rle>
#   [1]        1 1263583-1276694      *
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths



## Day 0
redpop_v2(
    atacF[1], chipF[1], sample_OCRs_Day0[41],
    atac.thres=30
)$res
# GRanges object with 0 ranges and 2 metadata columns:
#    seqnames    ranges strand | nearest_min_cov  max_atac
#       <Rle> <IRanges>  <Rle> |       <numeric> <numeric>
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Day 14
redpop_v2(
    atacF[2], chipF[2], sample_OCRs_Day14[46]
)$res
# GRanges object with 1 range and 2 metadata columns:
#      seqnames          ranges strand | nearest_min_cov  max_atac
#         <Rle>       <IRanges>  <Rle> |       <numeric> <numeric>
#  [1]        1 1273825-1274061      * |        -261.423    63.893
#  -------
#  seqinfo: 1 sequence from an unspecified genome; no seqlengths
