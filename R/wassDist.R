#' Compute the Wasserstein distance between to finite distributions.
#'
#' @param a,b Vectors representing probability distributions.
#' @param distMat Cost matrix.
#' @param p cost exponent.
#' @return The Wasserstein distance.
#' @export
wassDist <- function(a, b, distMat, p = 1){
  tr <- transport::transport(a, b, costm = distMat^p)
  tr <- as.matrix(tr)
  sum(distMat[tr[, 1:2]] * tr[, 3])^(1/p)
}
