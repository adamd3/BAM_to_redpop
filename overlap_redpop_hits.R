#!/usr/bin/env Rscript

library(GenomicRanges)
library(redpop)
library(dendextend)
library(ChIPseeker)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

## Function arguments:
## 1. chip_atac_combined.tsv = tsv file matching ATAC-Seq and ChIP-Seq samples
## 2. res_dir = directory containing redpop results
## 3. outdir = directory for output files

in_df = args[1]
res_dir = args[2]
outdir = args[3]


##------------------------------------------------------------------------------
## Read and process redpop results
##------------------------------------------------------------------------------
chip_atac_df <- read.table(in_df, header=TRUE)
fnames <- paste0(chip_atac_df[,3],"_",chip_atac_df[,4])
chromosomes <- 1:22

redpop_files <- gsub(" ","",file.path(paste0(
    apply(expand.grid(fnames, chromosomes), 1, paste0, collapse="_chr"),
    ".narrow.redpop.rds")))

redpop_files <- redpop_files[file.exists(redpop_files)]

## subset by chromosome
# redpop_files <- redpop_files[grep("chr1\\.",redpop_files)]

redpop_hits <- lapply(redpop_files, readRDS)

names(redpop_hits) <- gsub("\\.narrow\\.redpop\\.rds","",redpop_files)

# names(redpop_hits[[1]][[1]])
# [1] "res"         "h3k27ac.pb"  "atac.pb"     "sel.range"   "covs.min"
# [6] "peakiness"   "smooth_covs"
## - `res` contains the hit (GenomicRanges object)
## - `sel.range` contains the full open chromatin region that was scanned

redpop_res <- lapply(redpop_hits, function(x){
    res <- lapply(x, function(y){ y$res })
    do.call("c", res)  ## combine Granges objects
})

no_hits <- data.frame(
    sample = gsub("_chr1", "", names(redpop_res)),
    no_hits = unlist(lapply(redpop_res, length))
)
no_hits$sample <- factor(no_hits$sample, levels = no_hits$sample)

hits_plot <- ggplot(no_hits, aes(x = sample, y = no_hits)) +
    geom_bar(stat="identity") +
    # geom_density(alpha=.2, fill="#FF6666") +
    xlab("Sample") +
    ylab("Number of hits") +
    theme_bw() +
    theme(
        text = element_text(size = cc1*1.5),
        title = element_text(size = cc1*1.5),
        # legend.position = "bottom",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = cc1*1.5),
        legend.text = element_text(size = cc1*1.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        # panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
            colour = "black", size=cc1*1.5, angle=45, hjust=1
        ),
        axis.text.y = element_text(colour = "black", size=cc1*1.5),
        axis.title.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.title.y = element_text(
            colour = "black", size=cc1*1.5, vjust = 2.5
        )
    )


ggsave(
    hits_plot, file = paste0(outdir, '/', 'redpop_no_hits_narrow.png'),
    device = "png", units = "in",
    width = 16, height = 8, dpi = 300
)


## Union of all hits
unionAll <- Reduce(union, unlist(redpop_res))
length(unionAll)
# [1] 3828


## get indices of overlaps of each sample with the Union of all
intersectsUnion <- lapply(redpop_res, function(x) {
    na.omit(findOverlaps(unionAll, x, select=c("all")))
    ## returns indices of overlaps in unionAll
})

sapply(intersectsUnion, length)


## make a data frame of overlaps
unionDf <- data.frame(unionAll)
unionDf$id  <- paste0(unionDf$seqnames, "_", unionDf$start) ## unique ID for each hit

overlapsDf <- data.frame(matrix(
    nrow = length(unionDf$id),
    ncol = length(names(redpop_hits))
))
rownames(overlapsDf) <- unionDf$id
colnames(overlapsDf) <- names(redpop_hits)

overlapsDf[is.na(overlapsDf)] <- 0

## populate rows with `1` if there's an overlap
for (i in seq_along(overlapsDf)){
    olidx <- as.data.frame(intersectsUnion[[i]])[,1]
    overlapsDf[olidx, i] <- 1
}

saveRDS(unionAll, file = paste0(outdir,"/unionAll.rds"))
saveRDS(overlapsDf, file = paste0(outdir,"/overlapsDf.rds"))



##------------------------------------------------------------------------------
## Annotate the union set based on distance to nearest TSS
##------------------------------------------------------------------------------
ensemblGenes <- loadDb("/home/ad866/rds/hpc-work/CF_project/ATACSeq_analysis/ensemblGenes270820.db")

annotUnion <- annotatePeak(
    unionAll, tssRegion=c(-2000, 2000), TxDb = ensemblGenes,
    annoDb = "org.Hs.eg.db", assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c(
        "Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"
    ),
    addFlankGeneInfo = TRUE, flankDistance = 10000,
    sameStrand = FALSE,
    overlap = "TSS"
)
annotUnion <- as.GRanges(annotUnion)
annotUnion$id <- paste0(
    seqnames(annotUnion), "_", start(annotUnion)
)
annotUnion$type <- ifelse(
    annotUnion$distanceToTSS <= 2000, "Promoter", "Distal"
)

identical(annotUnion$id, rownames(overlapsDf))
# [1] TRUE

overlapsDf$type <- annotUnion$type

table(overlapsDf$type)
# Distal Promoter
#   1041     2787



##------------------------------------------------------------------------------
## Histogram of number of samples per hit
##------------------------------------------------------------------------------

hist_df <- data.frame(
    id = rownames(overlapsDf),
    type = overlapsDf$type,
    no_samples = rowSums(overlapsDf[1:ncol(overlapsDf)-1])
)

hist_df$type <- as.factor(hist_df$type)

median(subset(hist_df, type=="Promoter")$no_samples)
# [1] 5
median(subset(hist_df, type=="Distal")$no_samples)
# [1] 3

hist_df$type <- factor(hist_df$type, levels=c("Promoter","Distal"))

rpLensPlot <- ggplot(hist_df, aes(x = no_samples, fill=type)) +
    geom_histogram(
        colour = "black", binwidth = 1, position = 'dodge'
    ) +
    # geom_density(alpha=.2, fill="#FF6666") +
    xlab("Number of samples") +
    ylab("Number of hits") +
    theme_bw() +
    theme(
        text = element_text(size = cc1*1.5),
        title = element_text(size = cc1*1.5),
        # legend.position = "bottom",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = cc1*1.5),
        legend.text = element_text(size = cc1*1.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        # panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.text.y = element_text(colour = "black", size=cc1*1.5),
        axis.title.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.title.y = element_text(
            colour = "black", size=cc1*1.5, vjust = 2.5
        )
    )


ggsave(
    rpLensPlot, file = paste0(outdir, '/', 'redpop_hit_lengths.png'),
    device = "png", units = "in",
    width = 10, height = 8, dpi = 300
)

# all hits (not split by promoter/distal):
# png(
#     paste0(outdir, "/hits_by_no_samples.png"),
#     width = 900, height = 800, pointsize = 20
# )
# hist(
#     rowSums(overlapsDf[1:ncol(overlapsDf)-1]),
#     main = "Number of samples per redpop hit",
#     breaks = 1:ncol(overlapsDf),
#     ylab = "Number of hits", xlab = "Number of samples")
# dev.off()



##------------------------------------------------------------------------------
## Hierarchically cluster samples on presence/absence matrix:
##------------------------------------------------------------------------------
hclust_redpop <- hclust(dist(t(overlapsDf[1:ncol(overlapsDf)-1])))
dend_redpop <- as.dendrogram(hclust_redpop)

## colour factor for dendrogram

large_disc_pal <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colpal_large <- unlist(
    mapply(brewer.pal, large_disc_pal$maxcolors, rownames(large_disc_pal)))
colpal_large[5:8] <- colpal_large[60:63] ## replace to avoid colour clashes


# colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colpal <- vector("character", length(names(redpop_hits)))
colpal[grep("CF001", names(redpop_hits))] <- colpal_large[1]
colpal[grep("CF002", names(redpop_hits))] <- colpal_large[2]
colpal[grep("CF003", names(redpop_hits))] <- colpal_large[3]
colpal[grep("CF004", names(redpop_hits))] <- colpal_large[14]
colpal[grep("CF005", names(redpop_hits))] <- colpal_large[5]
colpal[grep("CF006", names(redpop_hits))] <- colpal_large[6]
colpal[grep("CF007", names(redpop_hits))] <- colpal_large[7]
colpal[grep("CF008", names(redpop_hits))] <- colpal_large[8]
colpal[grep("CF009", names(redpop_hits))] <- colpal_large[9]
colpal[grep("CF010", names(redpop_hits))] <- colpal_large[10]
colpal[grep("CF011", names(redpop_hits))] <- colpal_large[11]
colpal[grep("CF012", names(redpop_hits))] <- colpal_large[12]
colpal[grep("CF013", names(redpop_hits))] <- colpal_large[13]
colpal[grep("HV", names(redpop_hits))] <- colpal_large[14]

colpal2 <- colpal[order.dendrogram(dend_redpop)]
labels_colors(dend_redpop) <- colpal2
dend_redpop <- dend_redpop %>% set("branches_lwd", 4)


png(
    paste0(outdir, "/hclust_redpop.png"),
    width = 900, height = 800, pointsize = 20
)
par(mar = c(12,4,1,1))
plot(dend_redpop)
colored_bars(
    colpal,
    dend_redpop,
    rowLabels = "Individual"
)
dev.off()



##------------------------------------------------------------------------------
## Compare variation within individual vs variation between individuals
##------------------------------------------------------------------------------
individuals <- colnames(overlapsDf)[1:(ncol(overlapsDf)-1)]
individuals <- as.factor(sub("\\_.*", "", individuals))

table(individuals)
# CF001 CF002 CF003 CF004 CF005 CF006 CF007 CF008 CF009 CF010 CF011 CF012 CF013
#     1     2     2     2     3     1     3     3     2     3     3     2     3
# HV001 HV005 HV006 HV007
#     1     1     1     1

multisampled <- names(table(individuals)[table(individuals)>=2])

overlapsDf_multi <- overlapsDf[1:(ncol(overlapsDf)-1)][individuals %in% multisampled]

multi_ol <- lapply(multisampled, function(x){
    overlapsDf_multi[,grep(x,colnames(overlapsDf_multi))]
})
multi_ol2 <- lapply(multi_ol, function(x){
    x <- x[rowSums(x)>0,] ## remove union hits not in this individual
    rowSums(x)==ncol(x)
})
multi_ol3 <- lapply(multi_ol2, function(x){
    sum(x)/length(x)
})


## get fraction of shared hits within individual:
multi_ol_within <- lapply(seq_along(multisampled), function(x){
    individual_hits <- names(multi_ol2[[x]])
    individual_ol <- overlapsDf[1:(ncol(overlapsDf)-1)][
        rownames(overlapsDf[1:(ncol(overlapsDf)-1)]) %in% individual_hits,
    ]
    nsamps <- as.numeric(table(individuals[individuals==multisampled[x]])[
        table(individuals[individuals==multisampled[x]])>0
    ])
    individual_only <- individual_ol[,grepl(multisampled[x],colnames(individual_ol))]
    other_only <- individual_ol[,!grepl(multisampled[x],colnames(individual_ol))]
    res <- lapply(1:1000, function(z){
        ## randomly select 2 columns from the same individual
        random_individual <- individual_only[sample(ncol(individual_only),2)]
        # random_individual <- random_individual[rowSums(random_individual)>0,]
        # rowSums(random_individual)==ncol(random_individual)
        random_individual[,1]==random_individual[,2]
    })
    res2 <- lapply(seq_along(res), function(y){
        sum(res[[y]])/length(res[[y]])
    })
    mean(unlist(res2))
})



## get fraction of shared hits between individuals:
multi_ol_between <- lapply(seq_along(multisampled), function(x){
    individual_hits <- names(multi_ol2[[x]])
    individual_ol <- overlapsDf[1:(ncol(overlapsDf)-1)][
        rownames(overlapsDf[1:(ncol(overlapsDf)-1)]) %in% individual_hits,
    ]
    nsamps <- as.numeric(table(individuals[individuals==multisampled[x]])[
        table(individuals[individuals==multisampled[x]])>0
    ])
    individual_only <- individual_ol[,grepl(multisampled[x],colnames(individual_ol))]
    other_only <- individual_ol[,!grepl(multisampled[x],colnames(individual_ol))]
    res <- lapply(1:1000, function(z){
        ## randomly remove one column from the individual and replace with
        ## a randomly selected column from all other individuals
        random_individual <- individual_only[sample(ncol(individual_only),1)]
        random_other <- other_only[sample(1:ncol(other_only),1)]
        random_combined <- cbind(random_individual,random_other)
        # random_combined <- random_combined[rowSums(random_combined)>0,]
        rowSums(random_combined)==ncol(random_combined)
    })
    res2 <- lapply(seq_along(res), function(y){
        sum(res[[y]])/length(res[[y]])
    })
    mean(unlist(res2))
})


oldf <- data.frame(
    individual = multisampled,
    shared_within <- unlist(multi_ol_within),
    shared_between <- unlist(multi_ol_between)
)
colnames(oldf) <- c("individual", "shared_within", "shared_between")

oldf_melt <- reshape2::melt(oldf)
oldf_melt$individual <- as.factor(oldf_melt$individual)

median(subset(oldf_melt, variable=="shared_within")$value)
# [1] 0.4302961
median(subset(oldf_melt, variable=="shared_between")$value)
# [1] 0.4027003

olindplot <- ggplot(oldf_melt, aes(x = individual, y = value, fill=variable)) +
    geom_bar(position="dodge", stat="identity", colour = "black", width=0.8) +
    xlab("Individual") +
    ylab("Fraction of hits shared") +
    scale_fill_manual(
        "",
        labels = c("Within same individual", "Between different individuals"),
        values = colpal_large) +
    theme_bw() +
    theme(
        text = element_text(size = cc1*1.5),
        title = element_text(size = cc1*1.5),
        legend.position = "top",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = cc1*1.5),
        legend.text = element_text(size = cc1*1.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        # panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
            colour = "black", size=cc1*1.5,hjust=1,angle=45
        ),
        axis.text.y = element_text(colour = "black", size=cc1*1.5),
        axis.title.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.title.y = element_text(
            colour = "black", size=cc1*1.5, vjust = 2.5
        )
    )


ggsave(
    olindplot, file = paste0(outdir, '/', 'overlap_within_between.png'),
    device = "png", units = "in",
    width = 9, height = 8, dpi = 300
)




##------------------------------------------------------------------------------
## Plot the lengths of RedPop hits:
##------------------------------------------------------------------------------
library(ggplot2)

cc1=12

redPopResGR <- unlist(as(redpop_res, "GRangesList"))
redPopResGR <- GenomicRanges::reduce(redPopResGR)

redpop_hit_lens <- width(redPopResGR)
redpop_df <- data.frame(
    length = redpop_hit_lens
)

rpLensPlot <- ggplot(redpop_df, aes(x = length)) +
    geom_histogram(colour = "black", fill = "red", binwidth = 40) +
    geom_density(alpha=.2, fill="#FF6666") +
    xlab("Length (bp)") +
    ylab("Number of hits") +
    theme_bw() +
    theme(
        text = element_text(size = cc1*1.5),
        title = element_text(size = cc1*1.5),
        # legend.position = "bottom",
        panel.grid.major = element_blank(),
        legend.title = element_text(size = cc1*1.5),
        legend.text = element_text(size = cc1*1.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        # panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.text.y = element_text(colour = "black", size=cc1*1.5),
        axis.title.x = element_text(
            colour = "black", size=cc1*1.5
        ),
        axis.title.y = element_text(
            colour = "black", size=cc1*1.5, vjust = 2.5
        )
    )


ggsave(
    rpLensPlot, file = paste0(outdir, '/', 'redpop_hit_lengths.png'),
    device = "png", units = "in",
    width = 8, height = 8, dpi = 300
)


##------------------------------------------------------------------------------
## Plot coverage of redpop hits
##------------------------------------------------------------------------------
sel_ranges <- sapply(redpop_hits[[1]], `[`, 4)
sel_ranges <- unlist(as(sel_ranges, "GRangesList"))


## Choose a hit for comparison across samples:
dim(overlapsDf)
# [1] 3828   35

dim(subset(overlapsDf, CF002_0_chr1==0 & CF002_14_chr1==1))
# [1] 1113   35
## (almost 1/3rd of the union set is in Day 14 but not in Day 0)


head(subset(overlapsDf, CF002_0_chr1==0 & CF002_14_chr1==1))
#           CF001_0_chr1 CF002_0_chr1 CF002_14_chr1 CF003_0_chr1 CF003_14_chr1
# 1_906745             1            0             1            0             0
# 1_1030265            0            0             1            0             0
# 1_1116065            0            0             1            1             0
# 1_1231945            1            0             1            1             0
# 1_1273785            1            0             1            1             0
# 1_1407145            1            0             1            1             0

# Select 1273785

head(subset(overlapsDf, CF001_0_chr1==1 & HV005_0_chr1==0))
#           CF001_0_chr1 CF004_0_chr1 CF007_14_chr1 CF010_14_chr1 CF011_30_chr1
# 1_817265             1            0             0             0             1
# 1_906745             1            1             0             1             1
# 1_1064265            1            1             0             0             0
# 1_1471625            1            0             0             0             0
# 1_1780425            1            1             1             1             1
# 1_2151385            1            0             1             0             1
#           HV005_0_chr1
# 1_817265             0
# 1_906745             0
# 1_1064265            0
# 1_1471625            0
# 1_1780425            0
# 1_2151385            0

# Select 1471625

## get the index of the selected hit:
grtest <- GRanges(
    seqnames = "1", ranges=IRanges(start=1273785,width=1),strand = c('*'))

selected_hit <- as.data.frame(findOverlaps(sel_ranges,grtest))[,1]
selected_hit
# [1] 9


png(
    paste0(outdir, 'redpop_AURKAIP_CF002_Day_14.png'),
    width = 1600, height = 1600, res = 200, units = "px", pointsize = 12
)
plot(
    redpop_hits[[3]][[selected_hit]], collapseTranscripts = TRUE, gb = "hg38"#,
    # biomart_filters=list(ensembl_gene_id = geneID)
)
dev.off()

# png(
#     paste0(outdir, 'redpop_ATAD3B_CF001.png'),
#     width = 1600, height = 1600, res = 200, units = "px", pointsize = 12
# )
# plot(
#     redpop_hits[[3]][[selected_hit]], collapseTranscripts = TRUE, gb = "hg38"#,
#     # biomart_filters=list(ensembl_gene_id = geneID)
# )
# dev.off()
