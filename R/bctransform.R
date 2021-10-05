#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
#'
#'

bctransform = function(y, lambda=0) {
  if (lambda == 0L) {
    log(y)
  }
  else {
    (y^lambda - 1) / lambda
  }
}
