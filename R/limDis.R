#' Sample from the limiting distribution under the null.
#'
#' @param B number of samples to generate. Defaults to 500.
#' @param r vector of probabilities in the original problem.
#' @param distMat distance matrix in the original problem.
#' @param p cost exponent. Defaults to 1.
#' @return A vector of samples.
#' @export
limDisNull <- function(B = 500, r, distMat, p = 1){

  N <- nrow(distMat)

  # Matrix for constraints.
  A <- kronecker(rep(1, N), diag(1, N, N)) -
    kronecker(diag(1, N, N), rep(1, N))

  rhs <- as.vector(distMat)^p

  dir <- rep("<=", N^2)

  # Define the derivative. If Rcplex is available use it to solve the LP,
  # else use Rglpk.
  if(requireNamespace("Rcplex", quietly = TRUE)){
    deriv <- function(g){
      Rcplex::Rcplex(cvec = g, Amat = A, bvec = rhs, lb = -Inf, objsense = "max",
                     sense = "L")$obj
    }
  } else{
    deriv <- function(g){
      Rglpk::Rglpk_solve_LP(obj = g, mat = as.matrix(A),
                            bounds = list(
                              lower = list(ind = 1:N, val = rep(-Inf, N))),
                            dir = dir, rhs = rhs, max = T)$optimum
    }
  }

  # Covariance matrix for limiting G.
  Sigma <- - r %o% r + diag(r)

  # Create sample from limiting distribution.
  apply(MASS::mvrnorm(n = B, mu = rep(0, N), Sigma = Sigma), 1, deriv)^(1/p)
}

#' Sample from the limiting distribution under the null when the underlying
#' space is a grid.
#'
#' @param B Number of bootstrap samples to generate. Defaults to 500.
#' @param r vector of probabilities in the original problem. Is interpreted as
#'          a square matrix.
#' @param p cost exponent.
#' @return A vector of samples.
#' @export
limDisNullGrid <- function(B = 500, r, p=1){
  L <- sqrt(length(r))
  coord <- seq(from = 0, to = 1, length.out = L)
  distMat <- as.matrix(stats::dist(expand.grid(coord, coord)))
  limDisNull(B = B, r = r, distMat = distMat, p = p)
}



#' Sample from the limit distribution under the alternative.
#'
#' @param B Number of samples to generate.
#' @param r,s Number of counts giving the two samples.
#' @param distMat Distance matrix.
#' @param p Cost exponent. Defaults to 1.
#' @return A vector of samples.
#' @export
limDisAlt <- function(B = 1000, r, s, distMat, p = 1){

  N <- length(r)
  n <- sum(r)
  m <- sum(s)
  lambda <- m / (n + m)
  r <- r / n
  s <- s / m
  wd <- wassDist(a = r, b = s, distMat = distMat, p = p)

  # Define target function.
  foo <- function(z, h1, h2){
    if(!all(c(h1 - z * r, h2 - z * s) >= 0)){
      return(Inf)
    } else{
     (wassDist(a = h1 - z * r,
                   b = h2 - z * s, distMat = distMat, p = p)^p + z * wd^p)
    }
  }

  # Define the covariance functions.
  Sigma1 <- - r %o% r + diag(r)
  Sigma2 <- - s %o% s + diag(s)

  sample <- c()
  # Minimize the target function.
  for(k in 1:B){
    G <- MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = Sigma1) * (r > 0)
    H <- MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = Sigma2) * (s > 0)
    sample[k] <- stats::optimize(foo, lower = -1e8, upper = 1e8,
                          h1 = sqrt(lambda) * G, h2 = sqrt(1 - lambda) * H)$objective
  }

  return(wd^(1 - p) * sample / p)
}

#' m-out-of-n Bootstrap for the limiting distribution.
#'
#' @param r,s Vectors of counts giving the two samples.
#' @param distMat Distance matrix.
#' @param B The number of samples to generate. Defaults to 1000.
#' @param p Cost exponent. Defaults to 1.
#' @param gamma m = n^gamma. Defaults to 0.9.
#' @return A sample from the limiting distribution.
#' @export
limDisAltBoot <- function(r, s, distMat, B = 1000, p = 1, gamma = 0.9){
 n <- sum(r)
 m <- sum(s)
 nB <- n^gamma
 mB <- m^gamma

 ev <- function(){
   rB <- stats::rmultinom(n = 1, size = nB, prob = r / sum(r))
   sB <- stats::rmultinom(n = 1, size = mB, prob = s / sum(s))

   wassDist(rB / nB, sB / mB, distMat = distMat, p = p)
 }

 sample <- replicate(n = B, ev())

 sample <- sqrt(((nB * mB) / (nB + mB))) *
   (sample - wassDist(r / sum(r), s / sum(s), distMat, p = p))
}
