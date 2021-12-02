# BAM_to_redpop
Takes BAM files as input and returns redpop results.

## Introduction
redpop is a method developed by [Dr Ernest Turro](https://www.mountsinai.org/profiles/ernest-turro) for the identification of regulatory regions from a combination of ATAC-Seq and H3K27Ac ChIP-Seq libraries.

The method followed here was originally described in [Turro *et al.* (2020)](https://www.nature.com/articles/s41586-020-2434-2).

### Input files required
- De-duplicated ATAC-Seq BAM files.
- Normalised BigWig track files for both ATAC-Seq and ChIP-Seq.
- TSV file matching ATAC-Seq and ChIP-Seq samples.

### Software required
- [F-Seq](https://github.com/aboyle/F-seq)
- [bedtools](https://bedtools.readthedocs.io/)
- [R](https://www.r-project.org/)

### R packages required
- [redpop](https://gitlab.haem.cam.ac.uk/et341/redpop/)
- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [GenomeInfoDb](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html)
- [ChIPpeakAnno](https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)

## Overview of method
The process can be summarised as follows:
- generate BED files from de-duplicated ATAC-Seq BAMs
- call ATAC-Seq peaks using F-Seq
- identify open chromatin regions (OCRs) from ATAC-Seq peaks as follows:
    - resize all F-seq peaks calls to a minimum of 3.2 kb
    - merge overlapping peaks
    - truncate merged coordinates so that none extend beyond end of chromosome
- scan OCRs for regulatory elements with redpop

## Output
The below files will be produced.
- `BED` files, one per ATAC-Seq library
- `.fseqPeaks` - merged peak calls from [F-Seq](https://github.com/aboyle/F-seq) in [narrow peak format](https://software.broadinstitute.org/software/igv/node/270)
- `fseqOCRs.rds` contains the identified OCRs (list of [GRanges objects](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html), one per chromosome for each sample)
- `*_chr*.redpop.rds` files (one per chromosome for each pair of ATAC/ChIP samples) with redpop hits ([GRanges objects](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html))
