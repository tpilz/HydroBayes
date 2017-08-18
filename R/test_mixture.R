#' 1-D mixture distribution as test pdf
#' @param x \code{numeric}. Value for which the function will be evaluated.
#'
#' @return Density of the mixture distribution for candidate value \code{x}.
#'
#' @references Based on Matlab code from Appendix C of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
mixture <- function(x) {
  return(1/6*dnorm(x, -8, 1) + 5/6*dnorm(x, 10 ,1))
}
