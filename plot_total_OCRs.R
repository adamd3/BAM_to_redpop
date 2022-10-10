
## plot total OCRs per sample:
fseqOCRs_combined <- lapply(fseqOCRs, function(x) {
    unlist(as(x, "GRangesList"))
})

fseqOCRs_chr1 <- lapply(fseqOCRs, function(x) {
    unlist(as(x[[1]], "GRangesList"))
})

libIDs <- read.table(libs_file, header=TRUE)[c(3,4)]
libIDs <- paste0(libIDs[,1],"_",libIDs[,2])


OCR_counts_df <- data.frame(
    sample = libIDs,
    no_OCRs = as.numeric(lapply(fseqOCRs_chr1, length))
)

library(ggplot2)
cc1=12
hits_plot <- ggplot(OCR_counts_df, aes(x = sample, y = no_OCRs)) +
    geom_bar(stat="identity") +
    # geom_density(alpha=.2, fill="#FF6666") +
    xlab("Sample") +
    ylab("Number of OCRs") +
    ggtitle("Number of OCRs detected (chromosome 1)") +
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
    hits_plot, file = paste0(outdir, '/', 'fseq_no_OCRs.png'),
    device = "png", units = "in",
    width = 16, height = 8, dpi = 300
)
