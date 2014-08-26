#' Get affinity cutoff values that corresponds to a given P-value
#'
#' @title determine significance threshold for tRap affinities
#' @param pwm the position weight matrix
#' @param matrix.name name of the transfac matrix
#' @param pvalue.cutoff significance threshold for the TFBS hits
#' @return affinity cutoff that corresponds to the given P-value
#' @export
get.affinity.cutoff <- function(matrix.name, pwm, pvalue.cutoff) {
  l = ncol(pwm)
  params = exactgevparams[[l]][matrix.name,]
  return(exp(qgev(pvalue.cutoff, loc=params["location"], scale=params["scale"], shape=params["shape"], lower.tail = FALSE)))
}
