# This file contains several internally used functions

calc_ll <- function(res, lik, obs, abc_rho, abc_e, glue_shape) {
  if(lik == 1) {
    out <- log(res)
  } else if(lik == 2) {
    out <- res
  } else if(lik == 21) {
    m <- length(res)
    out <- - m/2 * log(2*pi) - sum(log(abc_e)) - 0.5 * sum( (get(abc_rho)(res, obs) / abc_e)^2 )
  } else if(lik == 22) {
    out <- min(abc_e - get(abc_rho)(res, obs)) # largest deviation of the distance value from tolerance counts
  } else if(lik == 31) {
    out <- glue_shape * log( max(1 - sum((res-obs)^2) / sum((obs-mean(obs))^2), 0) )
  } else stop("Argument 'lik' has to be one of {1,2,21,22}!")
  return(out)
}


prior_sample <- function(par.info, d, nc) {
  # generate prior sample
  if(par.info$initial == "uniform") {
    xt <- t(replicate(nc, runif(d, min = par.info$min, max = par.info$max)))
    if(d == 1) xt <- t(xt)
  } else if(par.info$initial == "normal") {
    xt <- mvrnorm(nc, mu = par.info$mu, Sigma = par.info$cov)
  } else if(par.info$initial == "latin") {
    lhs_sample <- randomLHS(nc, d)
    xt <- t(apply(lhs_sample, 1, function(x) qunif(x, min = par.info$min, max = par.info$max)))
    if(d == 1) xt <- t(xt)
  } else if(par.info$initial == "user") {
    xt <- par.info$val_ini
  } else {
    stop("Value 'initial' of argument list 'par.info' must be one of {'uniform', 'normal', 'latin', 'user'}!")
  }
  if(!is.matrix(xt))
    xt <- matrix(xt, ncol=d)
  if(!is.null(par.info$names))
    colnames(xt) <- par.info$names
  return(xt)
}

prior_pdf <- function(x, par.info, lik) {
  if (lik == 22) {
    lp <- 0 # set to zero if ABC method is used, TODO: not sure how to handle this?!
  } else {
    # calculate prior log-density
    if(par.info$prior == "uniform") {
      lp <- sum(dunif(x, min = par.info$min, max = par.info$max, log = T))
    } else if(par.info$prior == "normal") {
      lp <- sum(dmvnorm(x, mean = par.info$mu, sigma = par.info$cov, log = T))
    } else {
      lp <- get(par.info$prior)(x)
    }
  }
  return(lp)
}

check_outlier <- function(dens, x, N) {
  # mean log density of second half of chain samples as proxy for fitness of each chain
  proxy <- colMeans( dens )
  # calculate the Inter Quartile Range statistic (IQR method) of the chains
  quartiles <- quantile(proxy, probs = c(0.25,0.75))
  iqr <- abs(diff(quartiles))
  # identify outlier chains
  outliers <- which(proxy < quartiles[1] - 2*iqr)
  # outlier chains take state of one of the other chains (randomly sampled as in Vrugt, 2016 instead of best chain as in Vrugt et al., 2009)
  if(length(outliers) > 0) {
    new_states <- sample((1:N)[-outliers], length(outliers), replace = FALSE)
    x[outliers,] <- x[new_states,]
    dens[nrow(dens),outliers] <- dens[nrow(dens),new_states]
  }
  # output
  return(list(xt = x, p_x = dens[nrow(dens),], outliers = outliers))
} # EOF check_outlier

bound_par <- function(par, min, max, handle) {
  par_out <- par
  if(!is.null(min) && !is.null(max) && !is.null(handle)) {
    if(par < min || par > max) {
      if(handle == "bound") {
        par_out <- max(min(par, max), min)
      } else {
        stop("Value 'bound' of argument list 'par.info' must be one of {'bound'}!")
      }
    }
  }
  return(par_out)
}

calc_prop <- function(j, x, d, nc, delta, CR, nCR, pCR, c_val, c_star, p_g, beta0, par.info) {
  # initialise jump dx
  dx <- rep(0, d)

  ## sub-space of chain pairs
  # number of chain pairs to be used to calculate jump (equal selection probabilities)
  D <- sample(1:delta, 1, replace = TRUE)
  # sample chains for jump calculation: a != b != j
  samp <- sample((1:nc)[-j], D*2, replace = FALSE)
  a <- samp[1:D]
  b <- samp[(D+1):(2*D)]

  ## parameter sub-space
  # index of crossover value
  id <- sample(1:nCR, 1, replace = TRUE, prob = pCR)
  # d values from U[0,1]
  z <- runif(d)
  # derive subset A of selected dimensions (parameters)
  A <- which(z < CR[id])
  # how many dimensions/parameters sampled?
  d_star <- length(A)
  # A needs one value at least
  if(d_star == 0) {
    A <- which.min(z)
    d_star <- 1
  }

  # draw lambda values (as stated in text instead of Algorithm 5/6)
  lambda <- runif(d_star, min = -c_val, max = c_val)
  # jump rate
  gamma_d <- 2.38/sqrt(2*D*d_star)
  # select jump rate gamma: weighted random sample of gamma_d or 1 with probabilities 1-p_g and p_g, respectively
  g <- sample(x = c(gamma_d, 1), size = 1, replace = TRUE, prob = c(1-p_g, p_g))
  # small random disturbance
  zeta <- rnorm(d_star, sd=c_star)
  # jump differential
  jump_diff <- colSums(x[a,A, drop=F] - x[b,A, drop=F])
  # compute jump (differential evolution) for parameter subset
  dx[A] <- zeta + (1+lambda) * g * jump_diff
  # adjust jumping distance if desired
  dx[A] <- dx[A] * beta0

  # proposal
  xp <- x[j,] + dx

  # check and adjust proposal
  xp <- sapply(1:d, function(k) bound_par(xp[k], min = par.info$min[k], max = par.info$max[k], handle = par.info$bound))
  # make sure parameter names are retained
  if(!is.null(par.info$names)) names(xp) <- par.info$names

  return(list(xp=xp, id=id))
}

metropolis_acceptance <- function(p_xp, p_x, lik) {
  if(lik == 22) {
    # modified metropolis acceptance probability for ABC method, see Sadegh and Vrugt, 2014
    p_acc <- max(ifelse(p_xp >= p_x, 1, 0), ifelse(p_xp >= 0, 1, 0))
    accept <- p_acc == 1

  } else {
    # Metropolis acceptance probability
    p_acc <- min(max(-100, p_xp - p_x), 0)
    # accept or reject
    accept <- p_acc > log(runif(1))
  }
  return(accept)
}
