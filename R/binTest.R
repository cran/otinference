#' Two-sample test for multivariate data based on binning.
#'
#' @param x,y The two samples, rows are realizations.
#' @param L Number of bins in each dimension.
#' @param B Number of realizations of limiting distribution to simulate.
#' @return p-value.
#' @export
#' @examples
#' \dontrun{
#' x <- MASS::mvrnorm(n = 100, mean = c(0, 0), Sigma = diag(1, 2))
#' y <- MASS::mvrnorm(n = 100, mean = c(0, 0), Sigma = diag(2, 2))
#' pVal <- binWDTest(x, y)}
binWDTest <- function(x, y, L = 5, B = 100){
  breaks <- sm::binning(rbind(x, y), nbins = L)$breaks
  xbin <- sm::binning(x, breaks = breaks)
  ybin <- sm::binning(y, breaks = breaks)

  a <- as.vector(xbin$table.freq)
  b <- as.vector(ybin$table.freq)

  distm <- as.matrix(stats::dist(arrayInd(1:L^ncol(x), .dim = rep(L, ncol(x)))))

  wd <- wassDist(a / sum(a), b / sum(b), distMat = distm)

  sam <- limDisNull(B = B, r = (a / sum(a) + b / sum(b)) / 2, distMat = distm)

  pVal <- sum(sam >= sqrt(nrow(x) * nrow(y) / (nrow(x) + nrow(y))) * wd) / B
}
