#' Regulatory Element Detection using Patterns Of Peaks
#'
#' @param atac.bw Path to 40bp-binned, normalised ATAC-seq BigWig file, or a \code{GRanges} object
#' @param h3k27ac.bw Path to 40bp-binned, normalised H3K27ac ChIP-seq BigWig file, or a \code{GRanges} object
#' @param sel.range \code{GRanges} object with a single row specifying where to search for regulatory elements
#' @param smooth.cov.thres Covariance minima between ATAC and H3K27ac need to be below this threshold
#' @param atac.thres Minimum ATAC height
#' @param peakiness.flank Flank used to define whether, locally, ATAC > H3K27ac or vice versa
#' @param covariance.flank Flank used to compute local covariance between ATAC and H3K27ac
#' @param smoothing.flank Flank used to smooth the covariances before identifying minima
#' @param local.min.flank Flank used to assess whether a smoothed covariance is a local minimum
#' @param peakiness.net Include all segments within these many bp from covariance minimum that have the correct peakiness
#' @return A \code{redpop} object
#' @importFrom GenomicRanges setdiff GRanges reduce ranges mcols mcols<- seqnames GRangesList
#' @importFrom BiocGenerics start end width
#' @importFrom IRanges IRanges
#' @importFrom methods as
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb seqlevels<-
#' @export
redpop_v2 <- function(atac.bw, h3k27ac.bw, sel.range,
                  smooth.cov.thres=-1,
                  atac.thres=47,
                  peakiness.flank=400,
                  covariance.flank=400,
                  smoothing.flank=400,
                  local.min.flank=80,
                  peakiness.net=100) {
  stopifnot(length(sel.range)==1)

  if(class(atac.bw)!="GRanges")
    atac.gr=import(BigWigFile(atac.bw), selection = sel.range )
  if(class(h3k27ac.bw)!="GRanges")
    h3k27ac.gr=import(BigWigFile(h3k27ac.bw), selection = sel.range )


  ## AD: extend the scores from the (binned) bigwig across the range:
  atac.pb=rep(atac.gr$score,width(atac.gr))
  h3k27ac.pb=rep(h3k27ac.gr$score,width(h3k27ac.gr))

  smooth_covs=smooth_track(get_covariances(atac.pb, h3k27ac.pb, flank=covariance.flank), flank=smoothing.flank)
  names(atac.pb)=names(h3k27ac.pb)=names(smooth_covs)=start(sel.range):end(sel.range)

  covs.min=get_local_mins(smooth_covs, flank=local.min.flank)
  covs.min=covs.min & smooth_covs < smooth.cov.thres

  res=GRanges()
  pk=peakiness(atac.pb, h3k27ac.pb, peakiness.flank)
  mcols(res)=list(score=logical())
  if(sum(covs.min,na.rm=TRUE)>0) { # found at least one regulatory element, now finemap them:
                                   # find nearest region of ATAC peak vs H3K27ac trough ("peakiness")
    pos=as.numeric(names(which(covs.min)))
    covs.min.gr=GRanges(seqnames(sel.range), IRanges(start=pos,end=pos), min_cov=smooth_covs[which(covs.min)])

    pk.gr=GRanges(
        seqnames(sel.range), IRanges(start=start(sel.range):end(sel.range), end=start(sel.range):end(sel.range)), score=!is.na(pk) & pk>0
    )
    pk.gr.T=GenomicRanges::reduce(pk.gr[pk.gr$score])
    res=unique(unlist(as(nearest_or_within(pk.gr.T, covs.min.gr, peakiness.net), "GRangesList")))

    # record the minimum covariance of all local minima that are either the nearest or within peakiness.net bp of each stretch
    mcols(res)$nearest_min_cov=sapply(nearest_or_within(covs.min.gr, res), function(x) { min(x$min_cov) })
  }

  mcols(res)$max_atac=apply(as.data.frame(ranges(res)),1,function(x) { max(atac.pb[as.character(x[1]:x[2])]) })
 # res=subset(res,max_atac>atac.thres)
  res=res[res$max_atac>atac.thres,]
  # seqlevels(res)=paste("chr",c(1:22,"X"),sep="")

  structure(
    class="redpop",
      list(res=res, h3k27ac.pb=h3k27ac.pb, atac.pb=atac.pb, sel.range=sel.range, covs.min=covs.min, peakiness=pk, smooth_covs=smooth_covs))
}
