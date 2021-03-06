#' TEST Implementation of DREAM
#' @param prior A function(N,d) that draws N samples from a d-variate prior distribution.
#' Returns an N-by-d matrix.
#' @param pdf A function(prior) that calculates the log-density of the target distribution for given prior.
#' Returns an N-variate vector.
#' @param nc \code{numeric}. Number of chains evolved in parallel.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#' @param burnin \code{numeric}. Length of the burn-in period as portion of t (\code{burnin period = burnin * t}).
#' These samples from the Markov chain will not be included in the output. Default: 0.
#' @param adapt \code{numeric}. Length of the adaptation period as portion of t
#' (\code{adaptation period = adapt * t}). Will be used for the update of crossover probabilities and
#' the replacement of outlier chains. Default: 0.1.
#' @param updateInterval \code{integer}. Interval for crossover probability updates during the adaptation period.
#' @param delta \code{integer}. Maximum number of chain pairs used to generate the jump (default: 3).
#' @param c_val \code{numeric}. Lambda value is sampled from U[-c_val,c_val] (default: 0.1).
#' @param c_star \code{numeric}. Zeta value sampled from N[0,c_star]. Should be small compared to target
#' (i.e. in this case the normal) distribution. Default: 1e-12.
#' @param nCR \code{integer}. Length of vector with crossover probabilities for parameter subspace sampling (default: 3).
#' @param p_g \code{numeric}. Probability for gamma, the jump rate, being equal to 1. Default: 0.2.
#' @param beta0 \code{numeric}. Reduce jump distance, e.g. if the average acceptance rate is low (less than 15 \%).
#' \code{0 < beta0 <= 1}. Default: 1 (i.e. jump distance is not adjusted).
#' @param thin \code{integer}. Thinning to be applied to output in case of large \code{t}. See below.
#' @param checkConvergence \code{logical}. Shall convergence of the MCMC chain be checked? Currently implemented:
#' Calculating the Gelman-Rubin diagnostic. Takes a lot of time! Default: FALSE.
#' @param ncores \code{integer} specifying the number of cores to be used for parallel \code{pdf} evaluation
#' using the \code{\link[doMC]{doMC}} package (Linux only!). A value > 1 is only useful for complex \code{pdf}.
#' Default: 1.
#' @param verbose \code{logical}. Print progress bar to console? Default: TRUE.
#' @param DEBUG \code{logical}. Option enables further output for error and/or more in-depth analysis.
#' See below. Default: FALSE.
#'
#' @return \code{list} with named elements:
#'
#' \emph{chain}: a (1-burnin)*t/thin-by-d-by-nc array of parameter realisations for each iteration and Markov chain;
#'
#' \emph{density}: a (1-burnin)*t/thin-by-nc matrix of log-densities computed by \code{pdf} at each iteration for each Markov chain;
#'
#' \emph{runtime}: time of function execution in seconds;
#'
#' \emph{outlier}: a list with adapt*t vectors of outlier indices in nc (value of 0 means no outliers);
#'
#' \emph{AR}: a (1-burnin)*t/thin-by-nCR matrix giving the acceptance rate for each sample number and crossover value
#' (first element is NA due to computational reasons);
#'
#' \emph{CR}: a (1-burnin)*t/thin-by-nCR matrix giving the selection probability for each sample number and crossover value
#' (first element is NA due to computational reasons);
#'
#' \emph{R_stat}: if \code{checkConvergence == T} a (1-burnin)*t/thin-50-by-d matrix giving the Gelman-Rubin convergence diagnostic
#' (note that at least 50 observations are used to compute R_stat). Otherwise \code{NULL}.
#'
#' IF DEBUG == TRUE:
#'
#' \emph{DEBUG}: a list with the elements:
#'
#' \emph{J}: a t-by-nCR matrix of cumulated Euclidean jump distances during the Markov chain progressing;
#'
#' \emph{dx}: a t-by-nc-by-d array of jump proposals;
#'
#' \emph{dx_eff}: a t-by-nc-by-d array of accepted jumps;
#'
#' \emph{std}: a t-by-d matrix of standard deviations of the chain for each parameter.
#' Note: J_i = J_i-1 + sum( (dx_eff_i/std_i)^2 );
#'
#' \emph{gamma}: a t-by-nc matrix of jump rate values;
#'
#' \emph{lambda}: a t-by-nc-by-d array of lambda values;
#'
#' \emph{zeta}: a t-by-nc-by-d array of zeta values;
#'
#' \emph{jump_diff}: a t-by-nc-by-d array of jump differentials ( sum(X_a - X_b) ).
#'
#' @details To understand the notation (e.g. what is lambda, nCR etc.), have a look at Sect. 3.3
#' of the reference paper (see below).
#'
#' @references Code based on 'Algorithm 5' and 'Algorithm 6' of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @export
dream_test <- function(prior, pdf, nc, t, d,
                  burnin = 0, adapt = 0.1, updateInterval = 10, delta = 3, c_val = 0.1, c_star = 1e-12, nCR = 3, p_g = 0.2,
                  beta0 = 1, thin = 1, checkConvergence = FALSE, ncores = 1, verbose = TRUE, DEBUG = FALSE) {

  ### Argument checks ###
  if(nc <= delta*2)
    stop("Argument 'nc' must be > 'delta' * 2!")


  ### Initialisations ###
  # register cores
  registerDoMC(ncores)

  # track processing time
  timea <- Sys.time()

  # allocate chains (respecting thin) and density (all samples for outlier calculation) for output
  out_x <- array(NaN, dim=c((1-burnin)*t/thin,d,nc))
  p_x <- array(NaN, dim=c(t,nc))

  # Variables for crossover probability selection and acceptance monitoring
  J <- n_id <- n_acc <- rep(0, nCR)

  # vector of sampled crossover indices
  id <- rep(NA, nc)

  # crossover values and selection probabilities
  CR <- c(1:nCR)/nCR
  pCR <- rep(1, nCR)/nCR

  # initialize chains by sampling from prior
  xt <- prior(nc,d)
  if(!is.matrix(xt)) # if d=1
    xt <- matrix(xt, ncol=d)
  p_x[1,] <- pdf(xt)

  xp <- array(NA, dim = c(nc, d))

  # auxiliary variables
  out_AR <- out_CR <- array(NA, dim = c((1-burnin)*t/thin, nCR))
  if(checkConvergence)
    out_rstat <- array(NA, dim = c((1-burnin)*t/thin-50, d))
  else
    out_rstat <- NULL
  outl <- list(NULL)
  ind_out <- 0

  if(burnin == 0 && 1 %% thin == 0) {
    out_x[1,,] <- t(xt)
    out_AR[1,] <- rep(0,3)
    out_CR[1,] <- pCR
    ind_out <- ind_out+1
  }

  # DEBUG
  if(DEBUG) {
    out_J <- array(NA, dim=c(t, nCR))
    out_J[1,] <- J
    out_dx <- array(NA, dim=c(t, nc, d))
    out_dx[1,,] <- 0
    out_dx_eff <- array(NA, dim=c(t, nc, d))
    out_dx_eff[1,,] <- 0
    out_std <- array(NA, dim=c(t, d))
    out_std[1,] <- 0
    out_gamma <- array(NA, dim=c(t,nc))
    out_gamma[1,] <- 0
    out_lambda <- array(NA, dim=c(t, nc, d))
    out_lambda[1,,] <- 0
    out_zeta <- array(NA, dim=c(t, nc, d))
    out_zeta[1,,] <- 0
    out_jumpdiff <- array(NA, dim=c(t, nc, d))
    out_jumpdiff[1,,] <- 0
  }


  ### helper functions for vectorisation ###

## function generating the jumping distance for a specific chain
  jump <- function(j, nc, x, d, CR, nCR, pCR, delta, p_g, beta0, DEBUG) {
    # initialise
    dx <- rep(0, d)

    # prepare parameter subspace sampling (DREAM-specific extensin of DE-MC)
    # select delta (equal selection probabilities)
    D <- sample(1:delta, 1, replace = TRUE)
    # extract a != b != j
    # a <- R[draw[1:D]]
    # b <- R[draw[(D+1):(D*2)]]
    samp <- sample((1:nc)[-j], D*2, replace = FALSE)
    a <- samp[1:D]
    b <- samp[(D+1):(2*D)]
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

    # output
    output <- list(dx = dx, id = id)
    if(DEBUG) {
      lambda_t <- rep(0, d)
      lambda_t[A] <- lambda
      zeta_t <- rep(0, d)
      zeta_t[A] <- zeta
      jump_diff_t <- rep(0, d)
      jump_diff_t[A] <- jump_diff
      output[["DEBUG"]] <- list(gamma = g, lambda = lambda_t, zeta = zeta_t, jump_diff = jump_diff_t)
    }
    return(output)
  } # EOF jump

## function for pdf evaluation and acceptance/rejection of proposal
  stepfun <- function(prop, x_last, p_last, std_x) {

    # calculate density at proposal
    p_prop <- pdf(prop)

    # probability of acceptance (Metropolis acceptance ratio)
    # p_acc <- min(1, p_prop/p_last)
    p_acc <- min(max(-100, p_prop - p_last), 0)
    if(p_acc > log(runif(1))) { # larger than sample point from U[0,1]?
      out_prop <- prop # accept candidate parameters
      out_p_prop <- p_prop # accept density accordingly
      acc <- 1 # accpected
    } else {
      # retain previous values
      out_prop <- x_last
      out_p_prop <- p_last
      acc <- 0 # rejected
    }

    # jump distance (euclidean distance of parameter moves)
    dJ <- sum(((out_prop - x_last)/std_x)^2)

    # output
    return(list(xt = out_prop,
                p_xt = out_p_prop,
                dJ = dJ,
                acc = acc))
  } # EOF stepfun

## outlier detection and correction (DREAM-specific)
  check_outlier <- function(dens, x, N) {
    # mean log density of second half of chain samples as proxy for fitness of each chain
    proxy <- colMeans( dens )
    # calculate the Inter Quartile Range statistic (IQR method) of the chains
    quartiles <- quantile(proxy, probs = c(0.25,0.75))
    iqr <- diff(quartiles)
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



  ### Algorithm ###

  # progress indicator
  if(verbose)
    pb <- txtProgressBar(min = 2, max = t, style = 3)

  # evolution of nc chains
  for(i in 2:t) {
    # next progress message
    if (verbose)
      setTxtProgressBar(pb, i)

    # generate jump
    jump_out <- lapply(1:nc, function(j) jump(j, nc, xt, d, CR, nCR, pCR, delta, p_g, beta0, DEBUG))
    dx <- t(sapply(jump_out, function(x) x$dx))
    if(d == 1)
      dx <- t(dx)
    id <- sapply(jump_out, function(x) x$id)

    if(DEBUG) {
      out_gamma[i,] <- sapply(jump_out, function(x) x$DEBUG$gamma)
      if (d==1) {
        out_lambda[i,,] <- t(t(sapply(jump_out, function(x) x$DEBUG$lambda)))
        out_zeta[i,,] <- t(t(sapply(jump_out, function(x) x$DEBUG$zeta)))
        out_jumpdiff[i,,] <- t(t(sapply(jump_out, function(x) x$DEBUG$jump_diff)))
      } else {
        out_lambda[i,,] <- t(sapply(jump_out, function(x) x$DEBUG$lambda))
        out_zeta[i,,] <- t(sapply(jump_out, function(x) x$DEBUG$zeta))
        out_jumpdiff[i,,] <- t(sapply(jump_out, function(x) x$DEBUG$jump_diff))
      }
    }

    # calculate proposal
    xp <- xt + dx

    if(DEBUG) xt_prev <- xt # DEBUG

    # std within a chain for each parameter do calculate euclidean jump distance
    std_x <- apply(xt, 2, sd)

    # function evaluation and proposal acceptance/rejection
    if (ncores > 1)
      step_out <- foreach(j=1:nc) %dopar% stepfun(xp[j,], xt[j,], p_x[i-1,j], std_x)
    else
      step_out <- lapply(1:nc, function(j) stepfun(xp[j,], xt[j,], p_x[i-1,j], std_x))
    # distribute calculated values to variables
    xt <- t(sapply(step_out, function(x) x$xt))
    if(d == 1)
      xt <- t(xt)
    p_x[i,] <- sapply(step_out, function(x) x$p_xt)
    dJ_t <- sapply(step_out, function(x) x$dJ)
    acc_t <- sapply(step_out, function(x) x$acc)
    for (h in 1:length(id)) {
      J[id[h]] <- J[id[h]] + dJ_t[h]
      n_id[id[h]] <- n_id[id[h]] + 1 # how many times crossover of id used?
      n_acc[id[h]] <- n_acc[id[h]] + acc_t[h] # how many times a proposal was accpected?
    }

    # store for output (respect burn-in period and possible output thinning)
    if( (i > (burnin*t)) && (i %% thin == 0) ) {
      ind_out <- ind_out + 1

      # chain states
      out_x[ind_out,,] <- t(xt)

      # AR and CR
      out_AR[ind_out,] <- n_acc/n_id
      out_CR[ind_out,] <- pCR

      # convergence diagnostic (needs as least 50 observations)
      if(checkConvergence == T && ind_out > 50)
        out_rstat[ind_out-50,] <- R_stat(out_x[1:ind_out,,, drop = F])
    }

    # during adaptation period
    if (i <= (adapt*t)) {
      if (i%%updateInterval == 0) {
        # update selection probability of crossover by jump distance following Vrugt, 2016 instead of Vrugt et al., 2009 (different results?!)
        # favours larger jumps over smaller ones to speed up convergence
        if(any(J > 0)) {
          pCR <- J/n_id
          pCR[which(is.nan(pCR))] <- 1/nCR # if a specific n_id is zero, i.e. was not yet used
          pCR <- pCR/sum(pCR)
        }
      }

      # check for outliers and correct them
      check_out <- check_outlier(p_x[ceiling(i/2):i, ], xt, nc)
      xt <- check_out$xt
      p_x[i,] <- check_out$p_x
      outl[[i]] <- check_out$outliers
    }

    # DEBUG
    if(DEBUG) {
      out_J[i,] <- J
      out_dx[i,,] <- dx
      out_dx_eff[i,,] <- xt - xt_prev
      out_std[i,] <- std_x
    }

  } #  end chain processing

  # close progress bar
  if(verbose)
    close(pb)

  # track processing time
  timeb <- Sys.time()

  # prepare output
  output <- list(chain = out_x,
                 density = p_x[seq(burnin*t+thin, t, by=thin),],
                 runtime = timeb - timea,
                 outlier = outl,
                 AR = out_AR,
                 CR = out_CR,
                 R_stat = out_rstat)

  if(DEBUG)
    output[["DEBUG"]] <- list(J = out_J,
                              dx = out_dx,
                              dx_eff = out_dx_eff,
                              std = out_std,
                              gamma = out_gamma,
                              lambda = out_lambda,
                              zeta = out_zeta,
                              jump_diff = out_jumpdiff)


  return(output)
} # EOF
