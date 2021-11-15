# BAM_to_redpop
Takes BAM files as input and returns redpop results.

## Introduction
redpop is a method developed by [Dr Ernest Turro](https://www.mountsinai.org/profiles/ernest-turro) for the identification of regulatory regions from a combination of ATAC-Seq and H3K27Ac ChIP-Seq libraries. The method followed here was originally described in Turro et al. (2020): https://www.nature.com/articles/s41586-020-2434-2

### Input files required
- De-duplicated ATAC-Seq BAM files.
- Normalised BigWig track files for both ATAC-Seq and ChIP-Seq.
- Txt file containing a list of the ATAC-Seq library names (must match the BAM files)

### Software required
- [F-Seq](https://github.com/aboyle/F-seq)
- [bedtools](https://bedtools.readthedocs.io/)
- [R](https://www.r-project.org/)
- [redpop](https://gitlab.haem.cam.ac.uk/et341/redpop/)

### R packages required
- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [GenomeInfoDb](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html)

## Overview of method
The process can be summarised as follows:
- generate BED files from de-duplicated ATAC-Seq BAMs
- call ATAC-Seq peaks using F-Seq
- identify open chromatin regions (OCRs) from ATAC-Seq peaks
- scan OCRs for regulatory elements with redpop
