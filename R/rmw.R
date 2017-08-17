#' Random Walk Metropolis algorithm
#' @param prior A function that draws N samples from a d-variate prior distribution.
#' @param pdf A function that calculates the density of the target distribution for given parameters.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#'
#' @return \code{list} with named elements 'chain': a t x d matrix of parameter realisations for each iteration
#' of the Markov chain; and 'density': a t-variate vector of densities computed by \code{pdf} at each iteration
#' of the Markov chain.
#'
#' @references Code from 'Algorithm 1' of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @import MASS
rwm <- function(prior, pdf, t, d) {

  # Covariance matrix of proposal distribution (identity matrix)
  C <- (2.38/sqrt(d))^2 * diag(d)

  # allocate chain and density
  x <- matrix(NaN, nrow=t, ncol=d)
  p_x <- rep(NaN, t)

  # initialize chain by sampling from prior
  x[1,1:d] <- prior(1,d)
  p_x[1] <- pdf(x[1,1:d])

  # evolution of chain
  for (i in 2:t) {
    # candidate point: previous point + sample from d-variate normal proposal distribution with zero means and covariance matrix C
    xp <- x[i-1, 1:d] + mvrnorm(n=d, mu=rep(0,d), Sigma = C)
    # density of candidate point
    p_xp <- pdf(xp)
    # probability of acceptance (Metropolis acceptance ratio)
    p_acc <- min(1, p_xp/p_x[i-1])
    if(p_acc > runif(1)) { # larger than sample point from U[0,1]?
      x[i,1:d] <- xp # accept candidate parameters
      p_x[i] <- p_xp # accept density accordingly
    } else {
      x[i,1:d] <- x[i-1,1:d] # retain previous parameter values
      p_x[i] <- p_x[i-1] # retain previous density accordingly
    }
  }

  return(list(chain=x, density=p_x))
} # EOF
