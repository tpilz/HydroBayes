#' Differential Evolution Markov Chain (DE-MC) algorithm
#' @param prior A function(N,d) that draws N samples from a d-variate prior distribution.
#' Returns an N-by-d matrix.
#' @param pdf A function(prior) that calculates the density of the target distribution for given prior.
#' Returns an N-variate vector.
#' @param nc \code{numeric}. Number of chains evolved in parallel.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#'
#' @return \code{list} with named elements 'chain': a t-by-d-by-nc array of parameter realisations for each iteration
#' and Markov chain; and 'density': a t-by-nc matrix of densities computed by \code{pdf} at each iteration for each
#' Markov chain.
#'
#' @references Code from 'Algorithm 4' of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @import MASS
#' @export
de_mc <- function(prior, pdf, nc, t, d) {

  # Default jump rate
  gamma_RWM <- 2.38/sqrt(2*d)

  # allocate chains and density
  x <- array(NaN, dim=c(t,d,nc))
  p_x <- array(NaN, dim=c(t,nc))

  # initialize chains by sampling from prior
  x[1,,] <- t(prior(nc,d))
  p_x[1,] <- pdf(t(x[1,,]))

  # R-Matrix: index of chains for DE
  R <- array(NaN, dim=c(nc, nc-1))
  nc_t <- 1:nc
  for(i in 1:nc) R[i, 1:(nc-1)] <- nc_t[which(nc_t != i)]

  # evolution of nc chains
  for(i in 2:t) {
    # permute [1, ..., nc-1] nc times (to draw parameters a and b randomly later on)
    draw <- replicate(nc, sample(c(1:(nc-1)),2))
    # select jump rate gamma: weighted random sample of gamma_RWM or 1 with probabilities 0.9 and 0.1, respectively
    g <- sample(x = c(gamma_RWM, 1), size = 1, replace = TRUE, prob = c(0.9, 0.1))

    # proposal and accept/reject for each chain
    for (j in 1:nc) {
      # extract a != b != j
      a <- R[j, draw[1,j]]
      b <- R[j, draw[2,j]]
      if(a == b || a == j || b == j)
        stop("During chain selection something unexpected happened (a==b or a or b equal to current chain index)!")

      # create proposal via differential evolution
      xp <- x[i-1,,j] + g * (x[i-1,,a] - x[i-1,,b]) + rnorm(d, sd=1e-6)
      # calculate density at proposal
      p_xp <- pdf(xp)

      # probability of acceptance (Metropolis acceptance ratio)
      p_acc <- min(1, p_xp/p_x[i-1,j])
      if(p_acc > runif(1)) { # larger than sample point from U[0,1]?
        x[i,,j] <- xp # accept candidate parameters
        p_x[i,j] <- p_xp # accept density accordingly
      } else {
        x[i,,j] <- x[i-1,,j] # retain previous parameter values
        p_x[i,j] <- p_x[i-1,j] # retain previous density accordingly
      }
    }
  }

  return(list(chain=x, density=p_x))
} # EOF
