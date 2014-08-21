#' Compute local affinity p-values
#'
#' Use an empirical distribution to calculate a p-value for a given affinity
#' value using a set of sequences that are passed in.
#' @title Calculate local p-values
#' @param affnt
#' @param pwm
#' @param seq
#' @param Rmax
#' @param lambda
#' @param pseudo.count
#' @param gc.content
#' @param window.size
#' @param window.offset
#' @return p-value
#' @export
local.paffinity <- function(affnt, pwm, seq, Rmax = NULL, lambda = 0.7, pseudo.count = 1, gc.content = 0.5, window.size=NULL, window.offset=1) {
  # if window size is not null we slide a window across the affinities and sum
  # if window offset can be used to avoid too large overlaps of the windows
  # compute all local affinities
  local.aff = affinity(pwm, seq, Rmax, lambda, pseudo.count, gc.content, slide=TRUE)
  if (!is.null(window.size)) {
    local.aff = sapply(seq(1, (length(local.aff) - window.size), by=window.offset), function(idx) {
      sum(local.aff[idx:(idx + window.size)])
    })
  }
  # make an empirical distribution from it (put inf to get a min of 1/n)
  distr = ecdf(c(Inf, local.aff))
  # return the p-value
  return(1 - distr(affnt))
}
