#!/usr/bin/env Rscript

library(GenomicRanges)
library(redpop)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 1. atacDatDir = Directory containing ATAC-Seq BigWig files
## 1. chipDatDir = Directory containing ChIP-Seq BigWig files

if (!length(args)==3) {
    stop("Must supply exactly 3 arguments")
} else {
    in_df = args[1]
    atacDatDir = args[2]
    chipDatDir = args[3]
}



##------------------------------------------------------------------------------
## read in data
##------------------------------------------------------------------------------
fseqOCRs <- readRDS("fseqOCRs.rds")

chip_atac_df <- read.table(in_df, header=TRUE)
bwExt <- ".40.bw" ## file extension for BigWigs

atacF <- paste0(atacDatDir, "/", chip_atac_df$ATAC_sample, bwExt)
chipF <- paste0(chipDatDir, "/", chip_atac_df$Chip_sample, bwExt)




##------------------------------------------------------------------------------
## Run RedPop
##------------------------------------------------------------------------------
lapply(1:nrow(chip_atac_df), function(x){
    fseqOCR_selected <- fseqOCRs[redpop_df_1[x,]$ATAC_sample][[1]]
    redres <- lapply(seq_along(fseqOCR_selected), function(y){
        tmp <- tryCatch(
            redpop(atacF[x], chipF[x], fseqOCR_selected[y]),
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
    redres <- redres[!is.na(redres)]
    sample <- paste0(chip_atac_df$Individual,"_",chip_atac_df$Time)[x]
    saveRDS(redres, file = paste0(sample,".redpop.rds"))
})
