#' Load a PWM matrix from transfac
#'
#' Jaspar file for all matrices should have the format as provided by
#' jaspar: http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/non_redundant/all_species/matrix_only/matrix_only.txt
#' @title read a transfac matrix
#' @param file transfac matrix / jaspar file for all matrices
#' @return read.transfac returns a list of matrices with four rows for A, C, G, T
#' read.transfac.dir returns a list of matrices read from a directory,
#' named with the filenames of the matrix files.
#' @export
#' @seealso \code{\link{affinity}}
#' @author Matthias Heinig
read.transfac <- function(file) {
  matrices = list()
  lines = readLines(file)
  id = NULL
  accession = NULL
  for (line in lines) {
    cmd <- substr(line, 1, 2)
    if (cmd == "AC") {
      # Save any existing count matrix because we are starting a new one
      if (exists("count.matrix")) {
        mat <- count.matrix[,seq_len(cols)]
        attr(mat, "id") = id
        attr(mat, "accession") = accession
        matrices[[id]] <- mat
      }
      # Create a new count.matrix and set the accession number
      count.matrix <- matrix(nrow=4, ncol=100)
      cols <- 0
      rownames(count.matrix) = c("A", "C", "G", "T")
      accession <- strsplit(line, "\\s+")[[1]][[2]]
    } else if (cmd == "ID") {
      id <- strsplit(line, "\\s+")[[1]][[2]]
    } else if (grepl("[0-9][0-9]", cmd)) {
      # numbered lines are part of the count matrix
      idx <- as.numeric(cmd)
      # In the rare case that there are more than 100 positions in the motif
      if (ncol(count.matrix) < idx) {
        count.matrix.tmp <- count.matrix
        count.matrix <- matrix(nrow=4, ncol=idx)
        count.matrix[,seq_len(ncol(count.matrix.tmp))] <- count.matrix.tmp
      }
      items <- strsplit(line, "\\s+")[[1]]
      count.matrix[,idx] <- as.numeric(items[2:5])
      cols <- max(cols, idx)
    }
  }
  if (exists("count.matrix")) {
    mat <- count.matrix[,seq_len(cols)]
    attr(mat, "id") = id
    attr(mat, "accession") = accession
    matrices[[id]] <- mat
  }
  return(matrices)
}

#' Load a whole directory of TRANSFAC matrices
#'
#' read.transfac.dir returns a list of matrices read from a directory,
#' named with the filenames of the matrix files.
#' @title read a directory of transfac matrices
#' @param dir directory of transfac matrices
#' @return read.transfac returns a list of matrices with four rows for A, C, G, T
#' @export
#' @seealso \code{\link{read.transfac}} to read a single matrix
#' @author Matthias Heinig
read.transfac.dir <- function(dir) {
  files = list.files(dir)
  matrices = lapply(files, function(f) read.transfac(file.path(dir, f)))
  names(matrices) = sapply(matrices, attr, "id")
  return(matrices)
}

#' Normalize a PWM such that the columns sum to one to use with the
#' \pkg{seqLogo} package.
#'
#' @title normalize a PWM
#' @param matrix position specific count matrix
#' @return normalized count matrix
#' @export
normalize.pwm <- function(matrix) {
  return(matrix / rep(apply(matrix, 2, sum), each=4))
}

#' Load a generic position-specific count matrix
#'
#' @title read a PSCM
#' @param file a file containing position specific count matrices
#' @return returns a list of matrices with four rows for A, C, G, T
#' @export
#' @seealso \code{\link{read.transfac.dir}}
#' @author Matthias Heinig
read.pscm <- function (file) {
  matrices = list()
  lines = readLines(file)
  count.matrix = numeric()
  read = F
  first = T
  for (line in lines) {
    if (substr(line, 1, 1) == ">") {
      if (!first) {
        attr(count.matrix, "id") = id
        attr(count.matrix, "accession") = ac
        matrices[[length(matrices) + 1]] = count.matrix
        count.matrix = numeric()
      } else {
        first = F
      }
      ac = strsplit(substr(line, 2, nchar(line)), "\\s")[[1]][1]
      match = regexpr("/name='\\S+'", line)
      id = substr(line, match, match + attr(match, "match.length"))
      id = strsplit(id, "'")[[1]][2]
    }
    else {
      items = strsplit(line, "\\s+")[[1]]
      count.matrix = cbind(count.matrix, as.numeric(items[1:4]))
      rownames(count.matrix) = c("A", "C", "G", "T")
    }
  }
  # append last matrix
  attr(count.matrix, "id") = id
  attr(count.matrix, "accession") = ac
  matrices[[length(matrices) + 1]] = count.matrix
  names(matrices) = sapply(matrices, attr, "id")
  return(matrices)
}

#' Load a PWM matrix from Jaspar
#'
#' Jaspar file for all matrices should have the format as provided by
#' jaspar: http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/non_redundant/all_species/matrix_only/matrix_only.txt
#' @title read a jaspar matrix
#' @param file jaspar file for all matrices
#' @return returns a list of matrices with four rows for A, C, G, T
#' @export
#' @seealso \code{\link{read.transfac.dir}}
#' @author Matthias Heinig
read.jaspar <- function (file) {
  matrices = list()
  lines = readLines(file)
  count.matrix = numeric()
  read = F
  first = T
  for (line in lines) {
    if (substr(line, 1, 1) == ">") {
      if (!first) {
        rownames(count.matrix) = c("A", "C", "G", "T")
        attr(count.matrix, "id") = id
        attr(count.matrix, "accession") = ac
        matrices[[length(matrices) + 1]] = count.matrix
        count.matrix = numeric()
      } else {
        first = F
      }
      id = strsplit(substr(line, 2, nchar(line)), "\\s")[[1]][1]
      ac = strsplit(substr(line, 2, nchar(line)), "\\s")[[1]][2]
    }
    else {
      line = gsub("]", "", gsub("[", "",
        gsub("[ACGT]", "", line), fixed=T), fixed=T)
      items = strsplit(line, "\\s+")[[1]][-1]
      count.matrix = rbind(count.matrix, as.numeric(items))
    }
  }
  # append last matrix
  attr(count.matrix, "id") = id
  attr(count.matrix, "accession") = ac
  matrices[[length(matrices) + 1]] = count.matrix
  names(matrices) = sapply(matrices, attr, "id")
  return(matrices)
}

