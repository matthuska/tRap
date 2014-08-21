#' Use generalized extreme value distributions fitted for each matrix to
#' compute p-values of trap affinities.
#'
#' given the affinity compute the p value for the matrix and the length len
#' if params is NULL, precomputed parameters will be loaded other wise params
#' are c(shape, scale, loc) given to the generalized extreme value distribution
#' @title compute p-values for affinity values
#' @param affinity trap affinity value
#' @param matrix.name name of the transfac matrix (ignored if params not NULL)
#' @param len length of the region used to compute the affinity
#' @param params if params is NULL (default), precomputed parameters
#'        will be loaded otherwise params are c(shape, scale, loc) given to
#'        the generalized extreme value distribution
#' @param no.exact do not use the parameters that were learned for sequences of the exact same length as \code{len}
#' @return  p-value of the affinity score(s)
#' @references Thomas Manke et.al.
#' @export
#' @examples
#' \dontrun{
#' matrices = read.transfac.dir(dir)
#' seq = read some fasta or get from BSgenome package
#' matrix.name = "AP1.."
#' af = affinity(matrices[[matrix.name]], seq)
#' paffinity(af, matrix.name, 200)}
paffinity <- function(affinity, matrix.name, len, params = NULL, no.exact = FALSE) {
  # this function makes use of data stored in R/sysdata.rda (gevparams, exact)
  verbose = options()$verbose
  # first log transform the affinity
  laff = log(affinity) # this is base e
  # log transform the length
  if (is.null(params)){
    # load the parameters fitted previously
    # if there are paramters exactly for the length, then use it
    if (len <= length(exactgevparams) && !is.null(exactgevparams[[len]]) && !no.exact && matrix.name %in% rownames(exactgevparams[[len]])) {
      if (verbose) cat("p-value from gev with exact length\n")
      params = exactgevparams[[len]][matrix.name, c("shape", "scale", "location")]
    } else if (matrix.name %in% rownames(gevparams)) {
      if (verbose) cat("p-value from gev with interpolated length\n")
      # otherwise extrapolate the paramters
      llen = log10(len) # here it is base 10, ask Thomas Manke about it..
      # the parameters are log length dependent
      x = matrix(as.numeric(gevparams[matrix.name,1:6]), nrow=2)
      params = c(1, llen) %*% x
      # check if all of the parameters are valid
      # scale must be positive
      if (!is.na(params[2]) && params[2] <= 0) {
        warning("length dependent scale is <= 0: set to NA")
        params[2] = NA
      }
    }
  }
  if (is.null(params) || any(is.na(params))) {
    warning(paste("parameters missing for matrix", matrix.name))
    return(rep(NA, length(affinity)))
  }
  p = pgev(laff, loc=params[[3]], scale=params[[2]], shape=params[[1]], lower.tail = FALSE)
  names(p) = NULL
  return(p)
}
