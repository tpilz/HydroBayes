#' Differential Evolution Adaptive Metropolis (DREAM) algorithm
#' @param prior A function(N,d) that draws N samples from a d-variate prior distribution.
#' Returns an N-by-d matrix.
#' @param pdf A function(prior) that calculates the density of the target distribution for given prior.
#' Returns an N-variate vector.
#' @param nc \code{numeric}. Number of chains evolved in parallel.
#' @param t \code{numeric}. Number of samples from the Markov chain.
#' @param d \code{numeric}. Number of parameters.
#' @param burnin \code{numeric}. Length of the burn-in period as portion of t (\code{burnin period = burnin * t}).
#' Must be < 1! Default: 0.1.
#' @param delta \code{integer}. Maximum number of chain pairs used to generate the jump (default: 3).
#' @param c_val \code{numeric}. Lambda value is sampled from U[-c_val,c_val] (default: 0.3).
#' @param c_star \code{numeric}. Zeta value sampled from N[0,c_star]. Should be small compared to target
#' (i.e. in this case the normal) distribution. Default: 1e-6.
#' @param nCR \code{integer}. Length of vector with crossover probabilities for parameter subspace sampling (default: 3).
#' @param p_g \code{numeric}. Probability for gamma, the jump rate, being equal to 1. Default: 0.2.
#' @param beta0 \code{numeric}. Reduce jump distance, e.g. if the average acceptance rate is low (less than 15 \%).
#' \code{0 < beta0 <= 1}. Default: 1 (i.e. jump distance is not adjusted).
#' @param thin \code{integer}. Thinning to be applied to output in case of large \code{t}. See below.
#'
#' @return \code{list} with named elements:
#'
#' \emph{chain}: a t/thin-by-d-by-nc array of parameter realisations for each iteration and Markov chain;
#'
#' \emph{density}: a t/thin-by-nc matrix of densities computed by \code{pdf} at each iteration for each Markov chain;
#'
#' \emph{runtime}: time of function execution in seconds;
#'
#' \emph{outlier}: a list with t vectors of outlier indices in nc (value of 0 means no outliers);
#'
#' \emph{AR}: a t/thin-by-nCR matrix giving the acceptance rate for each sample number and crossover value
#' (first element is NA due to computational reasons);
#'
#' \emph{CR}: a t/thin-by-nCR matrix giving the selection probability for each sample number and crossover value
#' (first element is NA due to computational reasons);
#'
#' \emph{R_stat}: a t/thin-by-d matrix giving the Gelman-Rubin convergence diagnostic (first two elements are NA
#' due to computational reasons).
#'
#' @details To understand the notation (e.g. what is lambda, nCR etc.), have a look at Sect. 3.3
#' of the reference paper (see below).
#'
#' @references Code from 'Algorithm 5' of:
#'
#' Vrugt, J. A.: "Markov chain Monte Carlo simulation using the DREAM software package:
#' Theory, concepts, and MATLAB implementation." Environmental Modelling & Software, 2016, 75, 273 -- 316,
#' \url{http://dx.doi.org/10.1016/j.envsoft.2015.08.013}.
#'
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'
#' @import MASS
#' @export
dream <- function(prior, pdf, nc, t, d,
                  burnin = 0.1, delta = 3, c_val = 0.1, c_star = 1e-6, nCR = 3, p_g = 0.2,
                  beta0 = 1, thin = 1) {

  ### Argument checks ###
  if(nc <= delta*2)
    stop("Argument 'nc' must be > 'delta' * 2!")


  ### Initialisations ###
  # track processing time
  timea <- Sys.time()

  # allocate chains (respecting thin) and density (all samples for outlier calculation) for output
  x <- array(NaN, dim=c(t/thin,d,nc))
  p_x <- array(NaN, dim=c(t,nc))

  # Variables for crossover probability selection and acceptance monitoring
  J <- n_id <- n_acc <- rep(0, nCR)

  # R-Matrix: index of chains for DE
  R <- array(NaN, dim=c(nc, nc-1))
  nc_t <- 1:nc
  for(i in 1:nc) R[i, 1:(nc-1)] <- nc_t[which(nc_t != i)]

  # crossover values and selection probabilities
  CR <- c(1:nCR)/nCR
  pCR <- rep(1, nCR)/nCR

  # initialize chains by sampling from prior
  xt <- prior(nc,d)
  if(!is.matrix(xt))
    matrix(xt, ncol=d)
  p_x[1,] <- pdf(xt)

  if(1 %% thin == 0)
    x[i/thin,,] <- t(xt)

  # auxiliary variables
  out_AR <- out_CR <- array(NA, dim = c(t/thin, nCR))
  out_rstat <- array(NA, dim = c(t/thin, d))
  outl <- list(NULL)



  ### Algorithm ###

  # evolution of nc chains
  for(i in 2:t) {
    # permute [1, ..., nc-1] nc times (to draw parameters a and b randomly later on)
    draw <- apply(array(runif((nc-1)*nc), dim = c(nc-1, nc)), 2, order)
    # set jump vectors to zero
    dx <- array(0, dim = c(nc, d))
    # draw lambda values
    lambda <- runif(nc, min = -c_val, max = c_val)
    # std for each dimension
    std_x <- apply(xt, 2, sd)

    # proposal and accept/reject for each chain
    for (j in 1:nc) {
      # prepare parameter subspace sampling (DREAM-specific extensin of DE-MC)
      # select delta (equal selection probabilities)
      D <- sample(1:delta, 1, replace = TRUE)
      # extract a != b != j
      a <- R[j, draw[1:D,j]]
      b <- R[j, draw[(D+1):(D*2),j]]
      if(any(a == b) || any(a == j) || any(b == j))
        stop("During chain selection something unexpected happened (some a equals some b or some a or b equals j)!")
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

      # jump rate
      gamma_d <- 2.38/sqrt(2*D*d_star)
      # select jump rate gamma: weighted random sample of gamma_d or 1 with probabilities 1-p_g and p_g, respectively
      g <- sample(x = c(gamma_d, 1), size = 1, replace = TRUE, prob = c(1-p_g, p_g))

      # compute jump (differential evolution) for parameter subset
      dx[j, A] <- rnorm(d_star, sd=c_star) + (1+lambda[j]) * g * apply(xt[a,A, drop=F] - xt[b,A, drop=F], 2, sum)
      # adjust jumping distance if desired
      dx[j, A] <- dx[j, A] * beta0

      # compute proposal
      xp <- xt[j,] + dx[j,]
      # calculate density at proposal
      p_xp <- pdf(xp)

      # probability of acceptance (Metropolis acceptance ratio)
      p_acc <- min(1, p_xp/p_x[i-1,j])
      if(p_acc > runif(1)) { # larger than sample point from U[0,1]?
        xt[j,] <- xp # accept candidate parameters
        p_x[i,j] <- p_xp # accept density accordingly
        n_acc[id] <- n_acc[id] + 1
      } else {
        # retain previous values
        p_x[i,j] <- p_x[i-1,j]
        # set jump back to zero for pCR
        dx[j,] <- 0
      }

      # update jump distance for id
      J[id] <- J[id] + sum((dx[j,]/std_x)^2)
      # how many times crossover of id used?
      n_id[id] <- n_id[id] + 1

    } # end proposal

    # store for output
    if(i %% thin == 0) {
      # chain states
      x[i/thin,,] <- t(xt)

      # AR and CR
      out_AR[i/thin,] <- n_acc/n_id
      out_CR[i/thin,] <- pCR

      # convergence diagnostic
      if(i/thin > 2)
        out_rstat[i/thin,] <- R_stat(x[1:(i/thin),,, drop = F])
    }

    # during burn-in period
    if (i <= (burnin*t) ) {
      # update selection probability of crossover by jump distance following Vrugt, 2016 instead of Vrugt et al., 2009 (different results?!)
      # favours larger jumps over smaller ones to speed up convergence
      if(any(J > 0)) {
        pCR <- J/n_id
        pCR[which(is.nan(pCR))] <- 1/nCR # if a specific n_id is zero, i.e. was not yet used
        pCR <- pCR/sum(pCR)
      }
    }

    # outlier detection and correction (DREAM-specific)
    # mean log density of second half of chain samples as proxy for fitness of each chain
    proxy <- apply( log( p_x[ceiling(i/2):i, 1:nc] ), 2, mean)
    # calculate the Inter Quartile Range statistic (IQR method) of the chains
    quartiles <- quantile(proxy, probs = c(0.25,0.75))
    iqr <- diff(quartiles)
    # identify outlier chains
    outliers <- which(proxy < quartiles[1] - 2*iqr)
    # outlier chains take state of one of the other chains (randomly sampled as in Vrugt, 2016 instead of best chain as in Vrugt et al., 2009)
    if(length(outliers) > 0) {
      new_states <- sample((1:nc)[-outliers], length(outliers), replace = FALSE)
      xt[outliers,] <- xt[new_states,]
      p_x[i,outliers] <- p_x[i,new_states]
      outl[[i]] <- outliers # keep track of outliers
    } else
      outl[[i]] <- 0 # keep track of outliers

  } #  end chain processing

  # track processing time
  timeb <- Sys.time()

  # prepare output
  output <- list(chain = x,
                 density = p_x[seq(thin, t, by=thin),],
                 runtime = timeb - timea,
                 outlier = outl,
                 AR = out_AR,
                 CR = out_CR,
                 R_stat = out_rstat)

  return(output)
} # EOF
