#!/usr/bin/env Rscript

library(GenomicRanges)
library(redpop)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 2. atacDatDir = Directory containing ATAC-Seq BigWig files
## 3. chipDatDir = Directory containing ChIP-Seq BigWig files
## 4. bwExt = File extension for BigWigs
## 5. chr_select = Selected chromosome to scan

if (!length(args)==4) {
    stop("Must supply exactly 4 arguments")
} else {
    in_df = args[1]
    atacDatDir = args[2]
    chipDatDir = args[3]
    bwExt = args[4]
    chr_select = args[5]
}

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
    sample_OCRs <- subset(sample_OCRs, names(sample_OCRs)==chr_select)
    atac_file <- atacF[x]
    chip_file <- chipF[x]
    ## idenitifier for the sample:
    sample <- paste0(chip_atac_df$Individual,"_",chip_atac_df$Time)[x]
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
        chr_res <- lapply(seq_along(chr_ocrs)[1:300], function(y){
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
