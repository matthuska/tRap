#' The TRAP method to predict transcription factor binding affinity to DNA
#' sequences (Roider et al 2007)
#'
#' @name tRap-package
#' @aliases tRap
#' @title compute transcription factor binding site affinity
#' @docType package
#' @references Roider, H. G.; Kanhere, A.; Manke, T. & Vingron, M. Predicting transcription factor affinities to DNA from a biophysical model. Bioinformatics, 2007, 23, 134-141
#' @useDynLib tRap
#' @keywords package
#' @import stringr ROCR evd Rcpp
#' @examples
#' pwm = matrix(c(5, 4, 3, 1, 10, 12, 5, 3, 3, 5, 3, 10), nrow=4)
#' seq = "ACTGACGTGTGCACACGATGCTAGCTG"
#' affinity(pwm, seq)
NULL
#' @name jaspar
#' @title The JASPAR transcription factor database in the format expected by tRap
#' @docType data
#' @format A list where each element is a single JASPAR matrix and is named by
#' the JASPAR accession number
#' @references Sandelin A, Alkema W, Engstrom P, Wasserman WW, Lenhard B.
#' JASPAR: an open-access database for eukaryotic transcription factor binding profiles.
#' Nucleic Acids Res. 2004 Jan 1;32(Database issue):D91-4.
#' @source http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt
#' @keywords datasets
NULL
#' @name gevparams
#' @title Precalculated generalized extreme value distribution parameters for
#' the JASPAR and TRANSFAC database matrices.
#' @docType data
#' @format A matrix containing the parameters (columns) for each pwm matrix (rows)
#' @keywords datasets
NULL
#' @name exactgevparams
#' @title Precalculated generalized extreme value distribution parameters for
#' the JASPAR and TRANSFAC database matrices. Calculated for a specific set of
#' sequence lengths (100, 200, 500, 600, 800, 1000, 2000 and 5000 bp)
#' @docType data
#' @format A list where each element is a matrix of parameters (columns) for
#' each pwm (rows). The list index is the length of sequence for which the
#' parameters were fit.
#' @keywords datasets
NULL
