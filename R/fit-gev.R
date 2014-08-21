#' Fit a gev distribution for a pwm matrix.
#'
#' \code{sequences} is a list of character vectors. Each of these
#' character vectors contains promoter sequences of the same length,
#' since the parameters of the generalized extreme value distribution for
#' the pwm are dependent on the sequence length. For each set of
#' sequences of the same length, the GEV parameters are fit. Finally
#' these parameters are used in a linear model dependent on the logarithm
#' base 10 of the length of the sequence.
#'
#' @title Fit a gev distribution for a pwm matrix.
#' @param pwm position specific count matrix with 4 rows: A, C, G, T
#' @param sequences the promoter sequences to fit the model
#' @param gc.content GC content to be passed to the \code{\link{affinity}} function
#' @param both.strands compute affinity for both strands (default: TRUE)
#' @return   An object of class GevFit. It contains two elements:
#'  \item{params}{the length dependent parameters given as the regression
#'  coefficients shape0, shape1, scale0, scale1, loc0 and loc1. These can
#'  be used to compute the gev parameters of a sequence of length l as
#'  follows: shape = shape0 + shape1 * log10(l) etc.}
#' @author Matthias Heinig <heinig@@molgen.mpg.de>
#' @export
#' @import evd
#' @examples
#' pwm = matrix(c(5, 4, 3, 1, 10, 12, 5, 3, 3, 5, 3, 10), nrow=4)
#' sequences = lapply(c(100, 200, 300), function(l) sapply(1:100,
#' function(x) paste(c("A", "C", "G", "T")[sample(4, l, replace=TRUE)],
#' collapse="")))
#'
#' fit.gev(pwm, sequences)
fit.gev <- function(pwm, sequences, gc.content=0.5, both.strands=TRUE) {
  # check if sequences in one block all have equal length
  region.size = sapply(sequences, function(block) {
    l = sapply(block, nchar)
    stopifnot(all(l == l[1]))
    return(l[1])
  })

  names(sequences) = as.character(region.size)
  
  params = data.frame()
  for (r in region.size) {
    cat("fit for region size", r, "\n")
    # compute affinities for all sequences
    affinities = affinity(pwm, sequences[[as.character(r)]], gc.content=gc.content, both.strands=both.strands)
    # log them
    laff = log(affinities)

    # take only the ones > -Inf
    laff = laff[laff > -Inf]
    
    # fit the gev
    # this process is error prone, so we do it in a trycatch block
    fit = NULL
    tryCatch({fit = fgev(laff)},
             error = function(e) {
               cat("error while fitting\n")
               print(e)
             })

    # store the fit
    if (!is.null(fit)) {
      params = rbind(params, data.frame(region.size=r, loc=fit$estimate[1], scale=fit$estimate[2], shape=fit$estimate[3], convergence=fit$convergence, stringsAsFactors=F))
    } else {
      params = rbind(params, data.frame(region.size=r, loc=NA, scale=NA, shape=NA, convergence="failed", stringsAsFactors=F))
    }
  }
  # fit the linear models for the parameters dependent on the log region size
  # for each parameter we estimate an intercept and a slope

  loc = lm(loc ~ log10(region.size), data=params)
  scale = lm(scale ~ log10(region.size), data=params)
  shape = lm(shape ~ log10(region.size), data=params)

  p.loc = anova(loc)[1,"Pr(>F)"]
  p.scale = anova(scale)[1,"Pr(>F)"]
  p.shape = anova(shape)[1,"Pr(>F)"]
  
  # the final results are the six model parameters and p-values for the
  # quality of the fits
  res = list(
    params=data.frame(shape0=coef(shape)[1], shape1=coef(shape)[2],
      scale0=coef(scale)[1], scale1=coef(scale)[2],
      loc0=coef(loc)[1], loc1=coef(loc)[2],
      p.shape=p.shape, p.scale=p.scale, p.loc=p.loc),
    gev.by.length=params)
  class(res) = "GevFit"

  rownames(res$params) = NULL
  rownames(res$gev.by.length) = NULL
  
  return(res)
}

