#' Load a PWM matrix from transfac
#'
#' Jaspar file for all matrices should have the format as provided by
#' jaspar: http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/non_redundant/all_species/matrix_only/matrix_only.txt
#' @title read a transfac matrix
#' @aliases read.transfac read.transfac.dir read.jaspar
#' @param file transfac matrix / jaspar file for all matrices
#' @param dir directory of transfac matrices
#' @return read.transfac returns a list of matrices with four rows for A, C, G, T
#' read.transfac.dir returns a list of matrices read from a directory,
#' named with the filenames of the matrix files.
#' @export
#' @seealso \code{\link{affinity}}
#' @author Matthias Heinig
read.transfac <- function(file) {
  matrices = list()
  lines = readLines(file)
  count.matrix = numeric()
  read = F
  first = T
  id = NULL
  accession = NULL
  for (line in lines) {
    if (substr(line, 1, 2) == "ID") {
      id = gsub("\\s", "", substr(line, 3, nchar(line)))
    }
    if (substr(line, 1, 2) == "AC") {
      if (!first) {
        attr(count.matrix, "id") = id
        attr(count.matrix, "accession") = accession
        matrices[[length(matrices) + 1]] = count.matrix
      } else {
        first = F
      }
      accession = gsub("\\s", "", substr(line, 3, nchar(line)))
    }
    if (read) {
      if (substr(line, 1, 2) == "XX") {
        read = F
      } else {
        items = strsplit(line, "\\s+")[[1]]
        count.matrix = cbind(count.matrix, as.numeric(items[2:5]))
        rownames(count.matrix) = c("A", "C", "G", "T")
      }
    }
    if (substr(line, 1, 2) == "P0") {
      read = T
    }
  }
  # append last matrix
  attr(count.matrix, "id") = id
  attr(count.matrix, "accession") = accession
  matrices[[length(matrices) + 1]] = count.matrix
  names(matrices) = sapply(matrices, attr, "id")
  return(matrices)
}

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

