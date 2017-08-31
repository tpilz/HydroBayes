#' Gelman-Rubin convergence diagnostic
#'
#' Calculates diagnostic values to determine convergence of a MCMC algorithm.
#'
#' @param x \code{numeric} array of [no. of samples]-by-[no. of parameters]-by-[no. of chains] dimensions
#'
#' @return A [no. of parameters]-dimensional vector of Gelman-Rubin diagnostic values.
#'
#' @details The Gelman-Rubin value should be <= 1.2 for all parameters to declare convergence.
#' Otherwise the number of samples should be increased.
#'
#' @references Algorithm based on Gelman and Rubin (1992) obtained from (Eqs. 32-35):
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
R_stat <- function(x) {

  # Initialisations
  n <- dim(x)[3] # no. of chains
  t <- dim(x)[1] # no. of samples per chain
  d <- dim(x)[2] # no. of parameters

  # helper function to be applied over each parameter
  W_r <- function(y) {
    # a mean value for each chain
    x_mean <- 2 / (t-2) * colSums(y)
    # calculate squared differences of values from mean for each chain and sum up everything
    return(sum(apply(y, 1, function(z) (z - x_mean)^2)))
  }

  # within-cain variance for each parameter
  W <- 2 / (n * (t-2)) * apply(x[ceiling(t/2):t,,, drop=F], 2, W_r)

  # between chain variance for each parameter
  B_r <- function(y) {
    # a mean value for each chain
    x_mean <- 2 / (t-2) * colSums(y)
    # mean of chain-means
    x_mm <- mean(x_mean)
    # variance of chain means
    return( 1 / (2 * (n-1)) * sum((x_mean - x_mm)^2) )
  }
  B_t <- apply(x[ceiling(t/2):t,,, drop=F], 2, B_r)

  # estimate variance of the target distribution
  sigma <- (t-2) / t * W + 2*B_t

  # calculate Gelman-Rubin diagnostic by variance comparison for each parameter
  out <- sqrt( (n+1)/n * sigma/W - (t-2)/(n*t) )

} # EOF
