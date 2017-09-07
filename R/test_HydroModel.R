#' Simple one parameter linear storage model for testing purposes
#' @param forcing \code{numeric} vector of forcing (precipitation) values.
#' @param param \code{numeric}. Parameter of the model.
#'
#' @return \code{numeric} vector of responses (discharge).
#'
#' @details output_i = param * storage_i with\cr
#' storage_i = (forcing_i * Dt + storage_i-1) / (1 + param * Dt) and\cr
#' Dt = 1 and storage = 1 as initial condition.
#'
#' @references Code transferred from Matlab code obtained
#' during the Hydrocourse held in Luxembourg, Spring 2016.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#' @export

HydroModel <- function(forcing,param){
  if(param < 0)
    stop("Argument 'param' has to be >= zero!")
  if(any(forcing < 0))
    stop("Values of argument 'forcing' need to be >= zero!")

  # Initials
  Dt <- 1
  S0 <- 1
  Y=rep(NaN, length(forcing))

  # Apply linear store equation
  for (i in 1:length(forcing)) {
    S_Dt <- (forcing[i] * Dt + S0) / (1+param*Dt)
    Y[i] <- param*S_Dt
    S0 <- S_Dt
  }

  return(Y)
}
