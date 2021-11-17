#!/usr/bin/env bash


##------------------------------------------------------------------------------
## Parameters
##------------------------------------------------------------------------------
source 'params.txt'


## Required files:
## ATAC_libraries.txt = txt file containing ATAC-Seq library names (one per line)
## chip_atac_combined.tsv = tab delimited file matching ChIP-Seq and ATAC-Seq samples


##------------------------------------------------------------------------------
## Generate ATAC-Seq BED files for Fseq
##------------------------------------------------------------------------------
for library in $(awk '{print $1}' ATAC_libraries.txt); do
    bamToBed -i $bam_dir/$library.bam > $bed_dir/$library.bed
done



##------------------------------------------------------------------------------
## Run peak calling with Fseq
##------------------------------------------------------------------------------
## Fseq can be retrieved from: https://github.com/aboyle/F-seq

## It has also been downloaded, compiled and stored in:
## /rds-d3/project/who1000/rds-who1000-cbrc/F-seq/

## NOTE that Fseq outputs a separate file for each chromosome; hence
## also need to merge these.

for library in $(awk '{print $1}' ATAC_libraries.txt); do
    mkdir -p $fseq_dir/$library
    /rds-d3/project/who1000/rds-who1000-cbrc/F-seq/dist~/fseq/bin/fseq \
        -o $fseq_dir/$library -of npf $bed_dir/$library.bed
    ## Merge chromosomes:
    cat ${fseq_dir}/${library}/*[[:digit:]].npf \
        ${fseq_dir}/${library}/*X.npf ${fseq_dir}/${library}/*Y.npf \
        > ${fseq_dir}/$library.fseqPeaks
    rm -rf ${fseq_dir}/${library}
done



##------------------------------------------------------------------------------
## Identify Open Chromatin Regions (OCRs) for redpop
##------------------------------------------------------------------------------
Rscript --vanilla call_OCRs.R ATAC_libraries.txt $fseq_dir hg38 \
    2> call_OCRs.log &



##------------------------------------------------------------------------------
## Run redpop
##------------------------------------------------------------------------------
Rscript --vanilla run_redpop.R chip_atac_combined.tsv $atacDatDir $chipDatDir \
    2> run_redpop.log &
