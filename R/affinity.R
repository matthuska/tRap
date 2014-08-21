#' The TRAP method to predict transcription factor binding affinity to DNA
#' sequences (Roider et al 2007). Default is to compute the affinity for
#' the whole region by summing over all sites. If site specific affinities
#' are of interest set \code{slide} to true.
#'
#' @title compute transcription factor binding site affinity
#' @param pwm position specific count matrix with 4 rows: A, C, G, T
#' @param seq the sequence
#' @param Rmax the maximal affinity times the activity of the factor by
#'        default set to \code{exp(0.584 * ncol(pwm) - 5.66)}
#' @param lambda scaling parameter for the Berg and von Hippel mismatch energy
#' @param pseudo.count pseudo count to add to the PWM
#' @param gc.content GC content for the background model
#' @param slide if true the affinities for each sequence position will
#'        be returned
#' @param both.strands scan on both strands (default is TRUE for DNA) or only on the forward
#'        strand (for instance for RNA)
#' @return if slide is false, the function returns to overall affinity, otherwise
#'  for each of the \code{nchar(seq) - ncol(pwm) + 1} sequence positions
#'  the score of the site is returned.
#'  If the sequence contains gap characters "-" then the positions will be
#'  removed and if \code{slide} is true the positions with gaps will be
#'  filled with NA. This is usefull if you want to compare polymorphisms
#'  in promoter sequences.
#' @references Roider, H. G.; Kanhere, A.; Manke, T. & Vingron, M. Predicting transcription factor affinities to DNA from a biophysical model. Bioinformatics, 2007, 23, 134-141
#' @export
#' @examples
#' pwm = matrix(c(5, 4, 3, 1, 10, 12, 5, 3, 3, 5, 3, 10), nrow=4)
#' seq = "ACTGACGTGTGCACACGATGCTAGCTG"
#' affinity(pwm, seq)
affinity <- function(pwm, seq, Rmax=NULL, lambda=0.7, pseudo.count=1, gc.content=0.5, slide=FALSE, both.strands=TRUE) {
  # rows of the pwm are labeled A, C, G, T and cols are the sequence positions
  # if slide is true, the affinities for all positions are returned

  if (is.list(seq)) {
    seq = unlist(seq)
  }
  
  # check the sequence to avoid seg faults
  lapply(seq, function(x) {
    if (is.null(x) || is.na(x) || mode(x) != "character") {
      stop("sequence must be a character string of length >= ncol(pwm)")
    }
  })

  # handle gaps in the sequence (can arise when comparing aligned sequences)
  if (slide) {
    gap.pos = str_locate_all(seq, fixed("-"))
    gap.pos = lapply(gap.pos, function(x) x[,1])
    # Save the sequence lengths so that the list of return values can
    # be set to the correct length after we remove gaps
    seq.len.orig = nchar(seq)
  }

  # remove gaps
  seq = gsub("-", "", seq)

  if (any(nchar(seq) < ncol(pwm))) {
    stop("sequence must be a character string of length >= ncol(pwm)")
  }

  # set Rmax to default  if it's NULL
  Rmax = ifelse(is.null(Rmax), exp(0.584 * ncol(pwm) - 5.66), Rmax)
  
  # add a pseudo count to the matrix
  pwm = pwm + pseudo.count

  at.content = 1 - gc.content

  pwm = apply(pwm, 2, function(p) {
    maxAT = max(p[c(1,4)])
    maxCG = max(p[c(2,3)])
    if(maxAT > maxCG){
      transformed = c(
        log(maxAT / p[1]) / lambda, # A
        log((maxAT / at.content) * (gc.content / p[2])) / lambda, # C
        log((maxAT / at.content) * (gc.content / p[3])) / lambda, # G
        log(maxAT / p[4]) / lambda) # T
    } else{
      transformed = c(
        log((maxCG / gc.content) * (at.content / p[1])) / lambda, # A
        log(maxCG/p[2]) / lambda, # C
        log(maxCG/p[3]) / lambda, # G
        log((maxCG/gc.content)*(at.content/p[4])) / lambda) # T
    }
    if(maxAT == maxCG){
      transformed = log(maxAT / p) / lambda
    }
    return(transformed)
  })

  # now go for the c-code
  if (slide) {
    res = R_affinity_multi(pwm, ncol(pwm), seq, nchar(seq), Rmax, lambda, both.strands)
    # when we compute the affinities for all positions we have to fill the gaps
    # with NA in order to get result in the same length as seq - ncol(pwm) + 1
    gapped_empty = lapply(seq.len.orig, function(x) numeric(length=x - ncol(pwm) + 1))
    gapped = mapply(function(x, y, z) {
      if (length(y) > 0) {
        x[y] = NA
        x[-y] = z
      } else {
        x = z
      }
      return(x)
    }, gapped_empty, gap.pos, res, SIMPLIFY=FALSE)
    if (length(gapped) == 1) {
      res = gapped[[1]]
    } else {
      res = gapped
    }
  } else {
    res = R_affinity_sum_multi(pwm, ncol(pwm), seq, nchar(seq), Rmax, lambda, both.strands)
  }
  return(res)
}
