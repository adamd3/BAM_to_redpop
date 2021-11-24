#!/usr/bin/env bash


##------------------------------------------------------------------------------
## Parameters
##------------------------------------------------------------------------------
source 'params.txt'


## Required files:
## chip_atac_combined.tsv = txt file containing ATAC-Seq library names (one per line)
## chip_atac_combined.tsv = tab delimited file matching ChIP-Seq and ATAC-Seq samples


##------------------------------------------------------------------------------
## Generate ATAC-Seq BED files for Fseq
##------------------------------------------------------------------------------
for library in $(awk 'FNR>1 {print $1}' chip_atac_combined.tsv); do
    bamToBed -i $bam_dir/$library.bam > $bed_dir/$library.bed
done



##------------------------------------------------------------------------------
## Run peak calling with Fseq
##------------------------------------------------------------------------------
## Fseq can be retrieved from: https://github.com/aboyle/F-seq

## It has also been downloaded, compiled and stored in:
## /rds-d3/project/who1000/rds-who1000-cbrc/F-seq/

for library in $(awk 'FNR>1 {print $1}' chip_atac_combined.tsv); do
    mkdir -p $fseq_dir/$library
    /rds-d3/project/who1000/rds-who1000-cbrc/F-seq/dist~/fseq/bin/fseq \
        -o $fseq_dir/$library -of npf $bed_dir/$library.bed
    # ## Merge chromosomes:
    # cat ${fseq_dir}/${library}/*[[:digit:]].npf \
    #     ${fseq_dir}/${library}/*X.npf ${fseq_dir}/${library}/*Y.npf \
    #     > ${fseq_dir}/$library.fseqPeaks
    # rm -rf ${fseq_dir}/${library}
done



##------------------------------------------------------------------------------
## Identify Open Chromatin Regions (OCRs) for redpop
##------------------------------------------------------------------------------
Rscript --vanilla call_OCRs.R chip_atac_combined.tsv $fseq_dir hg38 \
    2> call_OCRs.log &



##------------------------------------------------------------------------------
## Run redpop
##------------------------------------------------------------------------------
# Rscript --vanilla run_redpop.R chip_atac_combined.tsv $atacDatDir $chipDatDir \
#     $bwExt 2> run_redpop.log &


## split the input file for parallel processing
tail -n +2 ../chip_atac_combined.tsv | split -l 4 - split_
for file in split_*
do
    head -n 1 ../chip_atac_combined.tsv > tmp_file
    cat "$file" >> tmp_file
    mv -f tmp_file "$file"
done

## `all` = scan all chromosomes; 1 = scan only chromosome 1, 2 = only 2, etc.
parallel Rscript --vanilla run_redpop.R {} $atacDatDir $chipDatDir \
    $bwExt all ::: split_*  2> run_redpop.log
