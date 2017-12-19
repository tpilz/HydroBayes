#' Adaptive Metropolis algorithm (scaling factor)
#' @param prior A function(N,d) that draws N samples from a d-variate prior distribution.
#' Returns an N-by-d matrix.
#' @param pdf A function(prior) that calculates the density of the target distribution for given prior.
#' Returns an N-variate vector.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#'
#' @return \code{list} with named elements 'chain': a t-by-d matrix of parameter realisations for each iteration
#' of the Markov chain; and 'density': a t-variate vector of densities computed by \code{pdf} at each iteration
#' of the Markov chain.
#'
#' @references Code from 'Algorithm 3' of (note error in line before the last therein):
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
am_sd <- function(prior, pdf, t, d) {

  # scaling factor of the covariance matrix
  s_d <- 2.38^2 / d

  # Covariance matrix of proposal distribution (identity matrix)
  C <- s_d * diag(d)

  # allocate chain and density
  x <- matrix(NaN, nrow=t, ncol=d)
  p_x <- rep(NaN, t)

  # initialize chain by sampling from prior
  x[1,1:d] <- prior(1,d)
  p_x[1] <- pdf(x[1,1:d])

  # first sample accepted
  accept <- 1

  # evolution of chain
  for (i in 2:t) {

    # adaptation of the scaling factor
    if ( (i %% 25 == 0) && (i < t/2) ) { # at every 25th iteration during burn-in (50% of t)
      # acceptance rate
      AR <- 100 * (accept / (i-1))
      # reduce scaling factor or ...
      if(AR<20) s_d <- 0.8 * s_d # 0.8 chosen arbitrarily, should be adapted (e.g. by trial and error)
      # ... increase scaling factor to obtain desired acceptance rate between 20 and 30%
      if(AR>30) s_d <- 1.2 * s_d # 1.2 chosen arbitrarily, should be adapted (e.g. by trial and error)
    }
    # Covariance adaptation
    C <- (s_d + 1e-4) * diag(d)

    # candidate point: previous point + sample from d-variate normal proposal distribution with zero means and covariance matrix C
    xp <- x[i-1, 1:d] + mvrnorm(mu=rep(0,d), Sigma = C)
    # density of candidate point
    p_xp <- pdf(xp)
    # probability of acceptance (Metropolis acceptance ratio)
    p_acc <- min(1, p_xp/p_x[i-1])
    if(p_acc > runif(1)) { # larger than sample point from U[0,1]?
      x[i,1:d] <- xp # accept candidate parameters
      p_x[i] <- p_xp # accept density accordingly
      accept <- accept +1 # increase count
    } else {
      x[i,1:d] <- x[i-1,1:d] # retain previous parameter values
      p_x[i] <- p_x[i-1] # retain previous density accordingly
    }
  }

  return(list(chain=x, density=p_x))
} # EOF
