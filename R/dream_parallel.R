#' Differential Evolution Adaptive Metropolis (DREAM) algorithm
#'
#' Implementation of the DREAM algorithm, including variations (such as DREAM(zs) and DEAM(ABC))
#' depending on certain argument settings.
#'
#' @param fun \code{character}. Name of a function(x, ...) which is evaluated at x being a
#' d-dimensional parameter vector. Returns a likelihood, log-likelihood, simulation values or summary statistics
#' depending on argument \code{lik}.
#' @param ... Additional arguments for \code{fun}.
#' @param lik \code{integer}. Flag specifying a likelihood function to be used, see details.
#' @param par.info A \code{list} of parameter information.
#' @param initial \code{character}. Method for prior sampling. One of: uniform - sampling from a uniform distribution;
#' normal - a (multivariate) normal distribution; latin - latin hypercube sampling; user - value(s) given by the user.
#' @param min \code{numeric}. A d-dimensional vector of minimum values for each parameter to sample from if
#' \code{initial} is 'uniform' or 'latin'. Also defines bounding region for the proposal, see \code{bound}.
#' @param max \code{numeric}. A d-dimensional vector of maximum values for each parameter to sample from if
#' \code{initial} is 'uniform' or 'latin'. Also defines bounding region for the proposal, see \code{bound}.
#' @param mu \code{numeric}. d-dimensional vector of parameter means if \code{initial} is 'normal'.
#' @param cov \code{numeric}. d-by-d positive-definite symmetric matrix of parameter covariances if
#' \code{initial} is 'normal'.
#' @param val_ini \code{numeric} nc-by-d-dimensional matrix of prior values if \code{initial} is 'user'.
#' @param bound \code{character}. What to do if the proposal parameter is outside the defined min-max limits.
#' One of: bound - proposal is set to min/max value if it is smaller/larger than the defined limit.
#' @param names \code{character} vector of length d with names for the parameters. These can be used within \code{fun}
#' (in this case, parameter input x of \code{fun} is a named vector) and will appear in the output list element 'chain'.
#' @param prior \code{character} specifying the prior pdf. One of: 'uniform' (prior is assumed to be uniformly distributed between
#' \code{min} and \code{max}, i.e. non-informative); 'normal' (prior assumed to be normally distributed with \code{mu}
#' and \code{cov}); name of a user-defined function(x) with x being a d-dimensional parameter vector and returning
#' the log-density of the d-variate prior distribution at x.
#' @param nc \code{numeric}. Number of chains evolved in parallel.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#' @param burnin \code{numeric}. Length of the burn-in period as portion of t (\code{burnin period = burnin * t}).
#' These samples from the Markov chain will not be included in the output. Default: 0.
#' @param adapt \code{numeric}. Length of the adaptation period as portion of t
#' (\code{adaptation period = adapt * t}). Will be used for the update of crossover probabilities (NB: useful
#' if some parameter are correlated as they will get a higher chance of being updated jointly) and
#' the replacement of outlier chains (NB: this can enhance convergence to the target distribution). Default: 0.1.
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
#' @param obs \code{numeric} vector of observations to be compared with output of \code{fun}. Only needed for some
#' realisations of \code{lik} (see details).
#' @param abc_rho \code{character}. Name of an ABC distance function(sim, obs) calculating the distance between
#' simulated ('sim', output of \code{fun}) and observed (\code{obs}) diagnostic values. Needed for ABC methods of
#' \code{lik}.
#' @param abc_e \code{numeric} vector of length of \code{obs} specifying the ABC tolerance value (\code{lik = 22})
#' or representing the standard deviations for each summary statistics (e.g. streamflow signatures) returned
#' by \code{fun} (\code{lik = 21}).
#' @param glue_shape \code{numeric} scalar value used for GLUE-based informal likelihood functions (see details).
#' @param lik_fun \code{character} specifying the name of a user-defined likelihood function(sim, obs) with sim
#' being the output of \code{fun} and obs a vector of corresponding observations (argument \code{obs} above).
#' @param past_sample \code{logical}. Shall proposal generation by sampling from past states be applied
#' (i.e., DREAM(zs))? Default: \code{FALSE}.
#' @param m0 \code{integer}. The initial archive size (needed if \code{past_sample == TRUE}). Default: \code{NULL}
#' (if \code{past_sample == TRUE} and \code{is.null(m0)} it will be set to \code{10*d}).
#' @param archive_update \code{integer}. Update the archive by appending the current state of each Markov chain
#' at every [\code{archive_update}]th iteration. Default: \code{NULL} (if \code{past_sample == TRUE} and
#' \code{is.null(archive_update)} it will be set to 10).
#' @param psnooker \code{numeric} value giving the probability for a snooker (instead of the conventional parallel
#' direction) jump. Might be useful for complex target distributions to enhance the diversity of proposals.
#' In such a case, a recommended value is 0.1. Otherwise, leave it at zero to prevent snooker updates (the default).
#' NOTE: Can only be applied if \code{past_sample == TRUE}!
#' @param ncores \code{integer} specifying the number of CPU cores to be used. If > 1, packages \code{\link[doMC]{doMC}}
#' (Linux only!) and \code{\link[parallel]{parallel}} are needed. Values > 1 only useful if \code{fun} is very
#' complex and computational demanding, otherwise multiple thread handling will cause function slowdown! Default: 1.
#' @param checkConvergence \code{logical}. Shall convergence of the MCMC chain be checked? Currently implemented:
#' Calculating the Gelman-Rubin diagnostic. Takes a lot of time! Default: FALSE.
#' @param verbose \code{logical}. Print progress bar to console? Default: TRUE.
#'
#' @return \code{list} with named elements:
#'
#' \emph{chain}: a (1-burnin)*t/thin-by-d+3-by-nc array of parameter realisations and the log-pdfs of the prior ('lp'),
#' likelihood ('ll'), and posterior ('lpost') distribution for each iteration and Markov chain;
#'
#' \emph{fx}: a (1-burnin)*t/thin-by-nc-by-[length of \code{fun}'s output] array of raw output of \code{fun},
#' corresponds with parameter realisations in chain;
#'
#' \emph{runtime}: time of function execution;
#'
#' \emph{outlier}: a list with adapt*t vectors of outlier indices in nc (value of 0 means no outliers);
#'
#' \emph{AR}: a (1-burnin)*t/thin-by-nCR matrix giving the acceptance rate for each sample number and crossover value
#' (first element is NA due to computational reasons);
#'
#' \emph{CR}: a (1-burnin)*t/thin-by-nCR matrix giving the selection probability for each sample number and crossover value
#' (first element is NA due to computational reasons).
#'
#' IF checkConvergence == TRUE:
#'
#' \emph{R_stat}: a (1-burnin)*t/thin-50+1-by-d matrix giving the Gelman-Rubin convergence diagnostic
#' (note that at least 50 observations are used to compute R_stat).
#'
#' @details To understand the notation (e.g. what is lambda, nCR etc.), have a look at Sect. 3.3
#' of the reference paper (see below).
#'
#' Likelihood options (argument \code{lik}):
#'
#' 1: \code{fun} returns the likelihood of parameter realisation x given some observations.
#'
#' 2: \code{fun} returns the log-likelihood.
#'
#' 21: ABC diagnostic model evaluation using a continuous fitness kernel, a variation of so-called noisy-ABC.
#'
#' 22: ABC diagnostic model evaluation using a Boxcar likelihood. \code{fun} has to return m diagnostic values
#' to be compared with observations. Requires arguments \code{obs}, \code{abc_rho}, and \code{abc_e}. Modification
#' in computation of Metropolis probability is used (Eq. 13 of Sadegh and Vrugt, 2014). Prior pdf set to zero!
#'
#' 31: GLUE with informal likelihood function based on the NSE: \code{glue_shape * log(NSE)}. \code{fun} needs to
#' return a time series of model simulations and \code{obs} should contain a time series of corresponding observations.
#' \code{glue_shape} affects the sampling: small values result in high parameter uncertainty, large values produce
#' narrow posterior pdfs of parameters. Should be in the range of 1-100 (experiment!).
#'
#' 99: A user-defined function, see argument \code{lik_fun}.
#'
#' @note If you want to use a non-implemented likelihood function with nuisance variables to be calibrated along
#' with the actual model parameters, it is suggested to implement the likelihood calculation directly into \code{fun}.
#'
#' @references Code based on:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' For ABC method see:
#'
#' Sadegh, M. and J. A. Vrugt: "Approximate Bayesian computation using Markov Chain Monte Carlo simulation:
#' DREAM_ABC". Water Resources Research, 2014, 50, 6767 -- 6787, \url{http://dx.doi.org/10.1002/2014WR015386}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @import MASS
#' @import doMC
#' @import parallel
#' @import lhs
#' @export
dream_parallel <- function(fun, ..., lik = NULL,
                  par.info = list(initial = NULL, min = NULL, max = NULL, mu = NULL, cov = NULL, val_ini = NULL,
                                  bound = NULL, names = NULL, prior = "uniform"),
                  nc, t, d,
                  burnin = 0, adapt = 0.1, updateInterval = 10, delta = 3, c_val = 0.1, c_star = 1e-12, nCR = 3,
                  p_g = 0.2, beta0 = 1, thin = 1, obs = NULL, abc_rho = NULL, abc_e = NULL, glue_shape = NULL, lik_fun = NULL,
                  past_sample = FALSE, m0 = NULL, archive_update = NULL, psnooker=0, ncores = 1, checkConvergence = FALSE, verbose = TRUE) {

  ### Argument checks ###
  if(!past_sample && (nc <= delta*2) )
    stop("Argument 'nc' must be > 'delta' * 2!")
  if(past_sample & is.null(m0))
    m0 <- 10*d
  if(!past_sample & is.null(m0))
    m0 <- nc
  if(past_sample & is.null(archive_update))
    archive_update <- 10
  if(ncores > 1)
    registerDoMC(ncores)

  # track processing time
  timea <- Sys.time()

  # Variables for crossover probability selection and acceptance monitoring
  J <- n_id <- n_acc <- rep(0, nCR)

  # vector of sampled crossover indices
  id <- rep(NA, nc)

  # crossover values and selection probabilities
  CR <- c(1:nCR)/nCR
  pCR <- rep(1, nCR)/nCR

  # number of lines for the diagnostic output: depends on burnin and thin whereas the prior will always occur
  if(burnin > 0 || thin > 1)
    out_t <- (1-burnin)*t/thin + 1
  else
    out_t <- t

  # initialise output variables
  chain <- array(NA, dim = c(out_t, d + 3, nc))
  parnames <- if(is.null(par.info$names)) paste0("par", 1:d) else par.info$names
  dimnames(chain) <- list(NULL, c(parnames, "lp", "ll", "lpost"), NULL)
  out_AR <- out_CR <- array(NA, dim = c(out_t, nCR))
  if(checkConvergence)
    out_rstat <- array(NA, dim = c(out_t-50, d))
  outl <- list(NULL)

  # initialise archive: sample from prior
  z <- prior_sample(par.info, d, m0)
  # get current state of each chain from the archive (last nc samples)
  xt <- z[(m0-nc+1):m0,, drop=F]
  # calculate prior log-density
  lp <- apply(xt, 1, prior_pdf, par.info = par.info, lik = lik)

  # evaluate fun for prior value
  if(ncores > 1) {
    res_fun <- mclapply(1:nrow(xt), function(i) get(fun)(xt[i,], ...), mc.cores = ncores)
    #res_fun <- mclapply(1:nrow(xt), function(i) get(fun)(xt[i,], init_dir, period, flood_thresh), mc.cores = 8)
    res_fun <- matrix(unlist(res_fun), nrow=nc, byrow = T) # nc-by-[length of fun output]
  } else {
    res_fun <- apply(xt, 1, function(i) get(fun)(i, ...))
    #res_fun <- apply(xt, 1, function(i) get(fun)(i, init_dir, period, flood_thresh))
    if(!is.matrix(res_fun)) res_fun <- t(res_fun)
    res_fun <- t(res_fun) # nc-by-[length of fun output]
  }

  # output variable to collect (raw) function output
  fx <- array(NA, dim = c(out_t, nc, ncol(res_fun)))

  # calculate log-likelihood
  ll <- apply(res_fun, 1, calc_ll, lik = lik, obs = obs, abc_rho = abc_rho, abc_e = abc_e,
              glue_shape = glue_shape, lik_fun = lik_fun)

  # calculate posterior log-density
  lpost <- array(NA, dim = c(t,nc)) # monitor all lpost values for outlier identification
  lpost[1,] <- lp + ll

  # write into output variables
  chain[1,1:d,] <- t(xt)
  chain[1,"lp",] <- lp
  chain[1,"ll",] <- ll
  chain[1,"lpost",] <- lpost[1,]
  out_AR[1,] <- rep(0,3)
  out_CR[1,] <- pCR
  fx[1,,] <- res_fun
  ind_out <- 1

  # progress indicator
  if(verbose)
    pb <- txtProgressBar(min = 2, max = t, style = 3)

  # evolution of nc chains
  for(i in 2:t) {

    # next progress message
    if (verbose)
      setTxtProgressBar(pb, i)
    # index to store output
    if( (i > (burnin*t)) && (i %% thin == 0) )
      ind_out <- ind_out + 1
    # other initialisations
    dx <- array(0, dim=c(nc, d))
    lp <- rep(0, nc)
    ll <- rep(0, nc)

    # calculate proposals
    res_t <- lapply(1:nc, function(j) calc_prop(j, xt, d, nc, delta, CR, nCR, pCR, c_val, c_star, p_g, beta0, par.info, past_sample, z, psnooker))
    xp <- sapply(res_t, function(x) x$xp)
    if(!is.matrix(xp)) xp <- t(xp)
    xp <- t(xp) # nc-by-d as xt
    id <- sapply(res_t, function(x) x$id)
    # calculate prior log-density
    lp_xp <- apply(xp, 1, prior_pdf, par.info, lik)

    # evaluate fun for proposals
    if(ncores > 1) {
      res_fun_t <- mclapply(1:nc, function(j) get(fun)(xp[j,], ...), mc.cores = 8)
      #res_fun_t <- mclapply(1:nc, function(j) get(fun)(xp[j,], init_dir, period, flood_thresh), mc.cores = 8)
      res_fun_t <- matrix(unlist(res_fun_t), nrow=nc, byrow = T) # nc-by-[length of fun output]
    } else {
      res_fun_t <- apply(xp, 1, get(fun), ...)
      #res_fun_t <- apply(xp[1:3,], 1, get(fun), init_dir, period, flood_thresh)
      if(!is.matrix(res_fun_t)) res_fun_t <- t(res_fun_t)
      res_fun_t <- t(res_fun_t) # nc-by-[length of fun output]
    }
    # calculate log-likelihood
    ll_xp <- apply(res_fun_t, 1, calc_ll, lik, obs, abc_rho, abc_e, glue_shape, lik_fun)
    # calculate posterior log-density
    lpost_xp <- lp_xp + ll_xp

    # Metropolis acceptance
    accept <- sapply(1:nc, function(j) metropolis_acceptance(lpost_xp[j], lpost[i-1,j], lik))
    # update auxiliary variables
    dx[accept,] <- xp[accept,] - xt[accept,] # jumped distance
    xt[accept,] <- xp[accept,] # accept candidate parameters
    lp[accept] <- lp_xp[accept] # prior density
    ll[accept] <- ll_xp[accept] # likelihood
    lpost[i,accept] <- lpost_xp[accept] # accept density accordingly
    lpost[i,!accept] <- lpost[i-1,!accept] # keep former lpost value for rejected chains
    res_fun[accept,] <- res_fun_t[accept,] # raw function output

    ids <- table(id) # count occurences of specific ids
    n_id[as.numeric(names(ids))] <- n_id[as.numeric(names(ids))] + ids # no. of times id was used
    id_acc <- table(id[accept]) # count occurences of specific accepted ids
    n_acc[as.numeric(names(id_acc))] <- n_acc[as.numeric(names(id_acc))] + id_acc # count acceptances per nCR
    std_x <- apply(xt, 2, sd) # sd among chains
    jump_dist <- apply(dx, 1, function(x) sum( (x/std_x)^2 )) # Euclidean jump distances
    for(j in 1:nc) J[id[j]] <-  J[id[j]] + jump_dist[j]  # monitoring of jump distances


    ## update archive
    if( !is.null(archive_update) && (i %% archive_update == 0) )
      z <- rbind(xt, z)

    ## store for output (respect burn-in period and possible output thinning)
    if( (i > (burnin*t)) && (i %% thin == 0) ) {
      # chain states
      chain[ind_out,1:d,] <- t(xt)
      chain[ind_out,"lp",] <- lp
      chain[ind_out,"ll",] <- ll
      chain[ind_out,"lpost",] <- lpost[i,]

      # function output
      fx[ind_out,,] <- res_fun

      # AR and CR
      out_AR[ind_out,] <- n_acc/n_id
      out_CR[ind_out,] <- pCR

      # convergence diagnostic (needs as least 50 observations)
      if(checkConvergence == T && ind_out > 50)
        out_rstat[ind_out-50,] <- R_stat(out_x[1:ind_out,,, drop = F])
    }


    ## during adaptation period
    if (i <= (adapt*t)) {
      if (i%%updateInterval == 0) {
        # update selection probability of crossover by jump distance following Vrugt, 2016 instead of Vrugt et al., 2009 (different results?!)
        # favours larger jumps over smaller ones to speed up convergence
        if(any(J > 0)) {
          pCR <- J/n_id
          pCR[which(is.nan(pCR))] <- 1/nCR # if a specific n_id is zero, i.e. was not yet used
          pCR <- pCR/sum(pCR)
        }

        # check for outliers and correct them (not if DREAM_zs, i.e. archive sampling is activated)
        if(!past_sample) {
          check_out <- check_outlier(lpost[ceiling(i/2):i, ], xt, nc)
          xt <- check_out$xt
          lpost[i,] <- check_out$p_x
          outl[[i]] <- check_out$outliers
        }
      } # update interval
    } # adaptation period

  } #  end chain processing


  # close progress bar
  if(verbose)
    close(pb)
  # track processing time
  timeb <- Sys.time()

  # prepare output
  output <- list(chain = chain,
                 fx = fx,
                 runtime = timeb - timea,
                 outlier = outl,
                 AR = out_AR,
                 CR = out_CR)
  if(checkConvergence)
    output[["R_stat"]] <- out_rstat

  return(output)
} # EOF
