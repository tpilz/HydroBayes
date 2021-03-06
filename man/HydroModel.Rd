% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_HydroModel.R
\name{HydroModel}
\alias{HydroModel}
\title{Simple one parameter linear storage model for testing purposes}
\usage{
HydroModel(forcing, param)
}
\arguments{
\item{forcing}{\code{numeric} vector of forcing (precipitation) values.}

\item{param}{\code{numeric}. Parameter of the model.}
}
\value{
\code{numeric} vector of responses (discharge).
}
\description{
Simple one parameter linear storage model for testing purposes
}
\details{
output_i = param * storage_i with\cr
storage_i = (forcing_i * Dt + storage_i-1) / (1 + param * Dt) and\cr
Dt = 1 and storage = 1 as initial condition.
}
\references{
Code transferred from Matlab code obtained
during the Hydrocourse held in Luxembourg, Spring 2016.
}
\author{
Tobias Pilz \email{tpilz@uni-potsdam.de}
}
