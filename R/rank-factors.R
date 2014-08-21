#' Rank transcription factors according to changes in binding
#' affinity by sequence polymorphisms like SNPs
#'
#' Unify all cases of comparing affinities for a set of pairs of sequence
#' possible design choices are
#' - computing affinities for the two alleles
#'   - globally in a region around the SNP
#'   - locally in a sliding window of length w (length of each PWM)
#' - computing p-values of the affinities
#'   - using gobally fitted GEV distributions
#'   - empirical p-values from the surounding sequence
#' @title Rank transcription factors according to changes in binding
#'        affinity by sequence polymorphisms like SNPs
#' @aliases rank.factors.for.pairs
#' @param matrices list of transfac matrices
#' @param seq1 sequence allele 1
#' @param seq2 sequence allele 2
#' @param local compute affinities locally in a sliding window or over
#'    the whole sequence
#' @param gev use the fitted generalized extreme value distribution for
#'    p-values (if FALSE, empirical distribution over the rest of the
#'    sequence)
#' @param snp.pos position of the SNP in the sequences (should always be
#'    the same and only needed for local or empirical p-values (gev ==
#'    FALSE)
#' @param margins compute the affinity in a margin around the SNP for
#'    empirical p-values and local == FALSE
#' @param empirical.offset when using empirical p-values and local ==
#'    FALSE the offset can be used to determine the spacing of windows to
#'    measure the affinities over regions of same size as the margin
#' @param both.strands compute affinity to both strands
#' @return   A data.frame with the fields a1, a2, p1, p2, ratio, log.ratio, min.p,
#'  min.prod. Each row corresponds to a matrix, the a fields contain the
#'  affinities, the p fields the corresponding p-values of the matrix and
#'  the two sequences. Further the ratio of the p-values, the log10 ratio,
#'  the min and the product of the p-values are given. If local is true,
#'  the min is over all positions.
#' @export
#' @references Roider, H. G.; Kanhere, A.; Manke, T. & Vingron, M. Predicting transcription factor affinities to DNA from a biophysical model. Bioinformatics, 2007, 23, 134-141
#' @examples
#' data(jaspar)
#' pwm = jaspar[[1]]
#' #pwm = matrix(c(5, 4, 3, 1, 10, 12, 5, 3, 3, 5, 3, 10), nrow=4)
#' seq1 = "actgacgtgtgcaCacgatgctagctg"
#' seq2 = "actgacgtgtgcaTacgatgctagctg"
#' rank.factors(list(mypwm=pwm), seq1, seq2)
rank.factors <- function(matrices, seq1, seq2, local=F, gev=T, snp.pos=NULL, margins=NULL, empirical.offset=1, both.strands=TRUE) {
  if (local) { # we slide the matrix across the SNP
    # iterate over all matrices
    pair.res = t(sapply(matrices, function(mat) {
      # extract the sequence around the snp for allele 1
      s1 = substr(seq1, snp.pos - ncol(mat) + 1, snp.pos + ncol(mat) - 1)
      aff1 = affinity(mat, s1, gc.content=0.5, slide=TRUE, both.strands=both.strands)
      # extract the sequence around the snp for allele 2
      s2 = substr(seq2, snp.pos - ncol(mat) + 1, snp.pos + ncol(mat) - 1)
      aff2 = affinity(mat, s2, gc.content=0.5, slide=TRUE, both.strands=both.strands)
      # compute the pvalues
      if (gev) {
        p1 = paffinity(aff1, attr(mat, "id"), ncol(mat))
        p2 = paffinity(aff2, attr(mat, "id"), ncol(mat))
      } else {
        p1 = local.paffinity(aff1, mat, seq1, gc.content=0.5)
        p2 = local.paffinity(aff2, mat, seq2, gc.content=0.5)
      }
      # print(c(length(aff1), length(aff2), length(p1), length(p2)))

      # compare p-values for the two alleles
      ratio = p1 / p2
      max.diff = which.max(abs(log(ratio)))

      # instead of the ratio we try the probability logic where we simulate a
      # nand by taking the min or the product
      min.p = min(c(p1, p2))
      prod.p = p1 * p2
      min.prod = which.min(prod.p)
      return(c(a1=NA, a2=NA, p1=p1[max.diff], p2=p2[max.diff], ratio=ratio[max.diff], log.ratio=log10(ratio[max.diff]), min.p=min.p, min.prod=prod.p[min.prod]))
    }))
  } else { # we use the whole sequence provided to compute affinities
    # extract the sequence around the snp for both alleles
    s1 = seq1
    s2 = seq2
    window.size = NULL
    if (!gev || !(is.null(margins) || is.null(snp.pos))) {
      # if we do not use gev, we have longer sequences since we need
      # the flanking sequences for the p-values, but the affinity is
      # computed within the margins only
      stopifnot(!any(is.null(c(snp.pos, margins))))
      if (length(margins) == 1) { # if only one margin is given: symmetric
        margins = rep(margins, 2)
      }
      s1 = substr(seq1, snp.pos - margins[1], snp.pos + margins[2])
      s2 = substr(seq2, snp.pos - margins[1], snp.pos + margins[2])
      window.size = sum(c(1, margins))
    }
    # iterate over all matrices
    pair.res = t(sapply(matrices, function(mat) {
      # compute the affinities
      aff1 = affinity(mat, s1, gc.content=0.5, both.strands=both.strands)
      aff2 = affinity(mat, s2, gc.content=0.5, both.strands=both.strands)
      # compute the pvalues
      if (gev) {
        p1 = paffinity(aff1, attr(mat, "id"), nchar(s1))
        p2 = paffinity(aff2, attr(mat, "id"), nchar(s2))
      } else {
        p1 = local.paffinity(aff1, mat, seq1, gc.content=0.5, window.size=window.size, window.offset=empirical.offset)
        p2 = local.paffinity(aff2, mat, seq2, gc.content=0.5, window.size=window.size, window.offset=empirical.offset)
      }
      # print(c(aff1, aff2, p1, p2))

      # compare p-values for the two alleles
      ratio = p1 / p2

      # instead of the ratio we try the probability logic where we simulate a
      # nand by taking the min or the product
      min.p = min(c(p1, p2))
      prod.p = p1 * p2
      return(c(a1=aff1, a2=aff2, p1=p1, p2=p2, ratio=ratio, log.ratio=log10(ratio), min.p=min.p, min.prod=prod.p))
    }))
  }
  return(pair.res)
}

rank.factors.for.pairs <- function(matrices, sequences, pairs, ...) {
  # - matrices: list of transfac matrices
  # - sequences: list of sequences as generated by read.fasta named by the
  #   names used in the pairs matrix
  # - pairs: 2 column matrix with sequence names of related sequences (SNPs)
  #   if pairs == NULL it is assumed that there are only two sequences
  # - ... parameters passed on to rank.factors

  results = data.frame()
  # iterate over the pairs
  dummy = apply(pairs, 1, function(pair) {
    # get the sequences
    fa.entry1 = sequences[[pair[1]]]
    fa.entry2 = sequences[[pair[2]]]
    cat(paste("pair:", fa.entry1$desc, fa.entry2$desc, "\n", sep="\n"))
    pair.res = rank.factors(matrices, fa.entry1$seq, fa.entry2$seq, ...)
    results <<- rbind(results, data.frame(seq1=rep(pair[1], nrow(pair.res)), seq2=rep(pair[2], nrow(pair.res)), matrix=rownames(pair.res), pair.res))
  })
  results[,"seq1"] = as.character(results[,"seq1"])
  results[,"seq2"] = as.character(results[,"seq2"])
  results[,"matrix"] = as.character(results[,"matrix"])
  return(results)
}

#' Evaluate Ranking
#'
#' Unknown
#' @title Evaluate Ranking
#' @param ranking
#' @param seq.name2matrix
#' @param roc.plot
#' @param add
#' @param col
#' @param predictors
#' @return Raw AUC value as well as ratio, min and product
evaluate.ranking <- function(ranking, seq.name2matrix, roc.plot=FALSE, add=FALSE, col="black", predictors=c("ratio", "min", "prod", "rawratio")) {
  known.matrix = sapply(strsplit(ranking[,"seq1"], "_"), function(items) items[2])
  predicted.matrix = ranking[,"matrix"]

  labels = apply(cbind(known.matrix, predicted.matrix), 1, function(pair) {
    pair[2] %in% seq.name2matrix[[pair[1]]]
  })

  # use the ratio
  pred.ratio = prediction(abs(ranking[,"log.ratio"]), labels) # large scores have positive class (TRUE)
  auc.ratio = performance(pred.ratio, "auc")@y.values[[1]]

  # use the min p
  pred.min = prediction(1 - ranking[,"min.p"], labels) # large scores have positive class (TRUE)
  auc.min = performance(pred.min, "auc")@y.values[[1]]

  # use the product
  pred.prod = prediction(1 - ranking[,"min.prod"], labels) # large scores have positive class (TRUE)
  auc.prod = performance(pred.prod, "auc")@y.values[[1]]

  # use the raw ratio
  if (!any(is.na(ranking[,c("a1", "a2")]))) {
    pred.raw = prediction(abs(log10(ranking[,"a1"] / ranking[,"a2"])), labels) # large scores have positive class (TRUE)
    auc.raw = performance(pred.raw, "auc")@y.values[[1]]
  } else {
    auc.raw = NA
  }
  added = F
  if (roc.plot) {
    if ("ratio" %in% predictors) {
      perf = performance(pred.ratio, "tpr", "fpr")
      plot(perf, main="ROC plot for SNPs in TFBS with known TF", add=add, lty="dashed", col=col)
      added = TRUE
    }
    if ("min" %in% predictors) {
      perf = performance(pred.min, "tpr", "fpr")
      plot(perf, main="ROC plot for SNPs in TFBS with known TF", add=added, lty="dotted", col=col)
      added = TRUE
    }
    if ("prod" %in% predictors) {
      perf = performance(pred.prod, "tpr", "fpr")
      plot(perf, main="ROC plot for SNPs in TFBS with known TF", add=added, lty="longdash", col=col)
      added = TRUE
    }
    if ("rawratio" %in% predictors) {
      perf = performance(pred.raw, "tpr", "fpr")
      plot(perf, main="ROC plot for SNPs in TFBS with known TF", add=added, lty="dotdash", col=col)
      added = TRUE
    }
  }
  return(c(auc.ratio=auc.ratio, auc.min=auc.min, auc.prod=auc.prod, auc.rawratio=auc.raw))
}
