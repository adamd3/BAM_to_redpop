#' Fill in gaps in a \code{GRanges} object with intervals with the specified score
#'
#' @param gr1 A \code{GRanges} object, possibly with gaps, with a numeric metadata column called score
#' @param gr2 A second \code{GRanges} object specifying where to fill in the gaps
#' @param score The default numeric value of score for the filled in ranges
#' @return A version of \code{gr1} in which gaps have been filled in with ranges having the default score
#' @importFrom GenomicRanges setdiff start width findOverlaps mcols nearest GRangesList
#' @importFrom stats cov
#' @export
#' @keywords internal
add_missing_ranges <- function(gr1,gr2, score=0) {
  sedi=setdiff(gr2, gr1)
  if(length(sedi)>0) {
    sedi$score=score
    return(sort(unlist(GRangesList(gr1,sedi))))
  } else {
    return(gr1)
  }
}

#' Extend \code{GRanges} object up to a specified minimum width
#'
#' @param gr A \code{GRanges} object
#' @param min_width Minimum number of bp
#' @return An extended \code{gr} annotated with maxmum/minimum max_atac/nearest_min_cov respectively in case of overlaps
#' @importFrom GenomicRanges mcols<- mcols width<- width start<- start
#' @importFrom S4Vectors queryHits
#' @export
#' @keywords internal
extend_until <- function(gr,min_width) {
  too.small=width(gr)<min_width
  start(gr[too.small])=start(gr[too.small])-min_width/2 + width(gr[too.small])/2
  width(gr[too.small])=min_width
  gr2=GenomicRanges::union(gr,gr)
  # annotate with highest ATAC and lowest nearest_min_cov; other annotations lost
  if(all(c("nearest_min_cov", "max_atac") %in% colnames(mcols(gr)))) {
    fo=findOverlaps(gr2,gr)
    spl=split(data.frame(fo), factor(queryHits(fo)))
    stopifnot(all(unique(queryHits(fo))==1:length(gr2)))
    mcols(gr2)$nearest_min_cov=rep(0,length(gr2))
    mcols(gr2)$max_atac=rep(0,length(gr2))
    mcols(gr2)[sapply(spl,nrow)==1,c("nearest_min_cov", "max_atac")] = mcols(gr)[do.call(rbind,spl[sapply(spl,nrow)==1])$subjectHits,c("nearest_min_cov","max_atac")]
    mcols(gr2)[sapply(spl,nrow)>1,"nearest_min_cov"] = sapply(spl[sapply(spl,nrow)>1], function(x) { min(mcols(gr)[x$subjectHits,"nearest_min_cov"]) })
    mcols(gr2)[sapply(spl,nrow)>1,"max_atac"] = sapply(spl[sapply(spl,nrow)>1], function(x) { max(mcols(gr)[x$subjectHits,"max_atac"]) })
  }
  return(gr2)
}

get_covariances <- function(vec1, vec2, flank) {
  c(rep(NA,flank),sapply((1+flank):(length(vec1)-flank), function(i) {
    cov(vec1[(i-flank):(i+flank)], vec2[(i-flank):(i+flank)])
  }),rep(NA,flank))
}

smooth_track <- function(vec, flank) {
  c(rep(NA,flank),sapply((1+flank):(length(vec)-flank), function(i) {
    mean(vec[(i-flank):(i+flank)],na.rm=TRUE)
  }), rep(NA,flank))
}

peakiness <- function(vec1, vec2, flank) { # where is vec1 peaky and vec2 troughy
  c(rep(NA,flank), sapply((1+flank):(length(vec1)-flank) , function(i) {
    vec1[i]/mean(vec1[(i-flank):(i+flank)])-vec2[i]/mean(vec2[(i-flank):(i+flank)])
  }), rep(NA,flank))
}

get_local_mins <- function(vec, flank) {
  c(rep(NA,flank), sapply((1+flank):(length(vec)-flank) , function(i) {
    vec[i] < min(vec[(i-flank):(i-1)]) && vec[i] < min(vec[(i+flank):(i+1)])
  }), rep(NA,flank))
}

#' Find nearest ranges and those within a certain distance
#'
#' @param gr1 A \code{GRanges} object
#' @param gr2 A \code{GRanges} object
#' @param net Number of bp
#' @return For each rangein gr2, the nearest range in gr1 and those within \code{net} bp
#' @importFrom GenomicRanges distance
#' @export
#' @keywords internal
nearest_or_within <- function(gr1, gr2, net=-1) { # for each range in gr1, return the nearest range in gr2 and all which are within net bp
  lapply(seq_along(gr2), function(i) {
    x=gr2[i]
    gr1[unique(sort(c(which(distance(x,gr1) < net), nearest(x,gr1))))]
  })
}
