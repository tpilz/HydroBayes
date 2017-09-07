# HydroBayes
Collection of R functions with sampling algorithms for use in a hydrological context employing Bayesian statistics

## DESCRIPTION
The package provides a number of functions with MCMC-based algorithms to sample probability distribution functions. The function are designed for the use in hydrological applications.

## WARNING
The package was devloped for personal use. Most of the sampling function contain the mere algorithm but are lacking any "cosmetics" such as argument checks etc. If you want to use the package, you should therefore take some care and you should know what you are doing.

I recommend the use of this package primarily to get acquainted with the algorithms.

## INSTALLATION

* command line installation (requiring the package devtools):

```R
install.packages("devtools") 
library(devtools)
install_github("tpilz/HydroBayes")
```

## FEEDBACK and BUGS
Feel free to comment via github issues: [>LINK<](https://github.com/tpilz/HydroBayes/issues).

Or clone the repository, make your changes and adjustments, and create a pull request.

## REFERENCES
Most of the code is based on the Matlab code from:

[Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 75, 273 – 316, doi: 10.1016/j.envsoft.2015.08.013, 2016.](http://dx.doi.org/10.1016/j.envsoft.2015.08.013)
