# This file contains several internally used functions

calc_ll <- function(res, lik, obs, abc_rho, abc_e, glue_shape, lik_fun) {
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
  } else if(lik == 99) {
    out <- get(lik_fun)(res, obs)
  } else stop("Argument 'lik' has to be one of {1,2,21,22}!")
  # avoid non-finite values, set to very small value
  out <- max(out, log(.Machine$double.xmin))
  # check output
  if(any(!is.finite(out))) stop("Calculated log-likelihood is not finite!")
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
    if(par.info$prior == "flat") { # non-informative prior, i.e. posterior ~ likelihood
      lp <- 0
    } else if(par.info$prior == "uniform") { # multiplicative prior, i.e. parameters must not be correlated!
      lp <- sum(dunif(x, min = par.info$min, max = par.info$max, log = T))
    } else if(par.info$prior == "normal") { # correlation among parameters respected via sigma
      lp <- dmvnorm(x, mean = par.info$mu, sigma = par.info$cov, log = T)
    } else {
      lp <- get(par.info$prior)(x)
    }
  }
  # avoid non-finite values, set to very small value
  lp <- max(lp, log(.Machine$double.xmin))
  # check output
  if(!is.finite(lp)) stop(paste0("Calculated log-prior is not finite: ", lp, "; x = ", x))
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
      } else if(handle == "reflect") {
        par_out <- par
        while(par_out < min || par_out > max) {
          if(par_out < min) par_out <- min + abs(diff(c(par_out, min)))
          if(par_out > max) par_out <- max - abs(diff(c(par_out, max)))
        }
      } else if(handle == "fold") {
        par_out <- par
        while(par_out < min || par_out > max) {
          if(par_out < min) par_out <- max - abs(diff(c(par_out, min)))
          if(par_out > max) par_out <- min + abs(diff(c(par_out, max)))
        }
      } else {
        stop("Value 'bound' of argument list 'par.info' must be one of {'bound', 'reflect', 'fold'}!")
      }
    }
  }
  return(par_out)
}

calc_prop <- function(j, x, d, nc, delta, CR, nCR, pCR, c_val, c_star, p_g, beta0, par.info, past_sample, z, psnooker) {

  # apply snooker or parallel direction update
  snooker <- FALSE
  if(runif(1) < psnooker) snooker <- TRUE

  ## sub-space of chain pairs
  # number of chain pairs to be used to calculate jump (equal selection probabilities)
  D <- sample(1:delta, 1, replace = TRUE)
  if(past_sample) {
    if(snooker) {
      # snooker update
      samp <- sample(1:nrow(z), 3, replace = FALSE)
      r1 <- z[samp[1],]
      r2 <- z[samp[2],]
      r3 <- z[samp[3],]
    } else {
      # parallel direction update
      # sample from state archive for jump calculation: a != b
      samp <- sample(1:nrow(z), D*2, replace = FALSE)
      a <- samp[1:D]
      b <- samp[(D+1):(2*D)]
      r1 <- z[a,, drop=F]
      r2 <- z[b,, drop=F]
    }
  } else {
    # sample from chains for jump calculation: a != b != j
    samp <- sample((1:nc)[-j], D*2, replace = FALSE)
    a <- samp[1:D]
    b <- samp[(D+1):(2*D)]
    r1 <- x[a,, drop=F]
    r2 <- x[b,, drop=F]
  }


  if(past_sample && snooker) {
    # index of crossover value
    # TODO: not needed if snooker but unsure what to do as id is needed later on (introduced error should be negligible?!)
    id <- sample(1:nCR, 1, replace = TRUE, prob = pCR)
    ## calc jump as snooker update
    # jump rate (according to ter Braak and Vrugt, 2008)
    g <- runif(1, min = 1.2, max = 2.2)
    # small random disturbance
    zeta <- rnorm(d, sd=c_star)
    # calculate projections of r2 and r3 to line  x - r1
    rp2 <- orth_proj(x[j,], r1, r2)
    rp3 <- orth_proj(x[j,], r1, r3)
    # jump
    dx <- zeta + g * (rp2-rp3)

  } else {

    ## parallel direction update
    # initialise jump dx
    dx <- rep(0, d)

    ## parameter sub-space
    # index of crossover value
    id <- sample(1:nCR, 1, replace = TRUE, prob = pCR)
    # d values from U[0,1]
    zt <- runif(d)
    # derive subset A of selected dimensions (parameters)
    A <- which(zt < CR[id])
    # how many dimensions/parameters sampled?
    d_star <- length(A)
    # A needs one value at least
    if(d_star == 0) {
      A <- which.min(zt)
      d_star <- 1
    }

    ## jump
    # draw lambda values (as stated in text instead of Algorithm 5/6)
    lambda <- runif(d_star, min = -c_val, max = c_val)
    # jump rate
    gamma_d <- 2.38/sqrt(2*D*d_star)
    # select jump rate gamma: weighted random sample of gamma_d or 1 with probabilities 1-p_g and p_g, respectively
    g <- sample(x = c(gamma_d, 1), size = 1, replace = TRUE, prob = c(1-p_g, p_g))
    # small random disturbance
    zeta <- rnorm(d_star, sd=c_star)
    # jump differential
    jump_diff <- colSums(r1[,A, drop=F] - r2[,A, drop=F])
    # compute jump (differential evolution) for parameter subset
    dx[A] <- zeta + (1+lambda) * g * jump_diff
    # adjust jumping distance if desired
    dx[A] <- dx[A] * beta0

  }


  # proposal
  xp <- x[j,] + dx

  # check and adjust proposal
  if(any(!is.finite(xp))) stop(paste0("Calculated proposal is not finite!", ifelse(snooker, " Snooker update was applied.", "")))
  if(!is.null(par.info$min) && !is.null(par.info$max)) {
    xp <- sapply(1:d, function(k) bound_par(xp[k], min = ifelse(length(par.info$min) > 1, par.info$min[k], par.info$min),
                                                   max = ifelse(length(par.info$max) > 1, par.info$max[k], par.info$max),
                                                   handle = par.info$bound))
  }
  # make sure parameter names are retained
  if(!is.null(par.info$names)) names(xp) <- par.info$names

  return(list(xp=xp, id=id))
}

metropolis_acceptance <- function(p_xp, p_x, lik, mt, p_zp) {
  if(lik == 22) {
    # modified metropolis acceptance probability for ABC method, see Sadegh and Vrugt, 2014
    p_acc <- max(ifelse(p_xp >= p_x, 1, 0), ifelse(p_xp >= 0, 1, 0))
    accept <- p_acc == 1

  } else if(mt>1) {
    # Modified Metropolis acceptance probability, Laloy and Vrugt, 2012
    p1 <- sum(exp(p_xp))
    p2 <- sum(exp(p_zp))
    if(p2 == 0) {
      accept <- FALSE
    } else {
    p_acc <- min(1, p1 / p2)
    accept <- p_acc > runif(1)
    }
  } else {
    # Metropolis acceptance probability
    p_acc <- min(max(-100, p_xp - p_x), 0)
    # accept or reject
    accept <- p_acc > log(runif(1))
  }
  if(!is.finite(accept)) stop(paste0("Calculate Metropolis acceptance is not finite: ", paste(accept, collapse = ", "),
                                     ";p_xp = ", paste(p_xp, collapse=", "),
                                     ", p_x = ", paste(p_x, collapse=", ")))
  return(accept)
}

# calculates the prthogonal projection of point p to the line going through points a and b
orth_proj <- function(a, b, p) {
  if (all(a == b)) {
    q <- p
  } else {
    # directional vector, i.e. the line between a and b
    u <- a-b
    # orthogonal projection (see analytical geometry basics)
    q  <- a + u * sum((p-a)*u) / sum(u*u)
  }
  # checkout output
  if(any(!is.finite(q))) stop(paste0("Result of orthogonal projection is not finite: ", paste(q, collapse = ", "),
                                     "; a = ", paste(a, collapse = ", "), ", b = ", paste(b, collapse = ", "), ", p = ", paste(p, collapse = ", ")))
  return(q)
}
