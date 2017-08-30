#' Multivariate t-/Student distribution with 60 deg. of freedom
#' @param x \code{numeric} vector with values of the parameters of the distribution or a or matrix
#' with rows corresponding to observations and columns to parameters defining the dimensionality of the problem.
#'
#' @return A vector of type \code{numeric} and length equal to the number of rows of \code{x}
#' containing the log of the t-/Student pdf at \code{x}.
#'
#' @references Based on Matlab code from Appendix C of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @import mvtnorm
#'
#' @export
t_distribution <- function(x) {

  # convert input to matrix if needed
  if(!is.matrix(x))
    x <- matrix(x, ncol=length(x))

  # dimensionality of target distribution
  d <- dim(x)[2]

  # degrees of freedom
  df <- 60

  # correlation matrix
  A <- 0.5 * diag(d) + array(0.5, dim = c(d,d))

  # calculate covariance matrix
  C <- array(NA, dim = c(d,d))
  for (i in 1:d) {
    for (j in 1:d) {
      C[i,j] = A[i,j] * sqrt(i*j)
    }
  }

  # compute log-density of multivariate Student distribution
  log_l <- dmvt(x, sigma = C, df = df, log = TRUE)

  return(log_l)
}
