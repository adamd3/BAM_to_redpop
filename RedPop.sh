#!/usr/bin/env bash


##------------------------------------------------------------------------------
## Parameters
##------------------------------------------------------------------------------
source 'params.txt'


## Required files:
## libraries.txt = txt file containing ATAC-Seq library names (one per line)


##------------------------------------------------------------------------------
## Generate BED files for Fseq
##------------------------------------------------------------------------------
for library in $(awk '{print $1}' libraries.txt); do
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

for library in $(awk '{print $1}' libraries.txt); do
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
## Identify Open Chromatin Regions (OCRs) for RedPop
##------------------------------------------------------------------------------

Rscript --slave --vanilla call_OCRs.R libraries.txt $fseq_dir hg38 > call_OCRs.Rout
