#' Identify TFBS exceeding the given P-value (or affinity) cutoff
#'
#' @title identify TFBS hits in the sequence(s)
#' @param seq character vector with sequences to scan
#' @param pwm the position weight matrix
#' @param matrix.name name of the transfac matrix
#' @param pvalue.cutoff significance threshold for the TFBS hits
#' @param affinity.cutoff score threshold (default NULL)
#' @return a matrix with 4 columns: seq, start, end, affinity containing the positions of significant hits
#' @export
scan.sequence <- function(seq, pwm, matrix.name, pvalue.cutoff, affinity.cutoff=NULL) {
  if (is.null(affinity.cutoff)) {
    affinity.cutoff = get.affinity.cutoff(matrix.name, pwm, pvalue.cutoff)
  }
  aff = affinity(pwm, seq, slide=T)
  if (length(seq) == 1) {
    aff = list(aff)
  }
  len = sapply(aff, length)
  clen = c(1, cumsum(len))
  aff = unlist(aff)
  hits = which(aff > affinity.cutoff)
  if (length(hits) > 0) {
    hits = cbind(seq=findInterval(hits, clen, all.inside=T, rightmost.closed=T), start=hits, end=hits + ncol(pwm), affinity=aff[hits])
  } else {
    hits = matrix(nrow=0, ncol=4)
    colnames(hits) = c("seq", "start", "end", "affinity")
  }
  return(hits)
}

