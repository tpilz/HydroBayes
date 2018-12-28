
# Script reproducing examples from Vrugt (2016) Environ. Modell. Softw.

library(MASS)
library(mvtnorm)
library(lhs)
library(HydroBayes)

setwd("/home/tobias/R/MyPackages/HydroBayes/examples")


### default values for dream
burnin = 0
adapt = 0.1
updateInterval = 10
delta = 3
c_val = 0.1
c_star = 1e-12
nCR = 3
p_g = 0.2
beta0 = 1
thin = 1
checkConvergence = FALSE
verbose = F
obs = NULL
abc_rho = NULL
abc_e = NULL
glue_shape = NULL
lik_fun = NULL
past_sample = T
m0 = NULL
archive_update = NULL
psnooker=0.1
ncores = 1
mt = 1

### example from Fig. 4 ###

# uniform prior between -10 and 10
prior <- function(N,d) t(replicate(N, runif(d, -10, 10)))

# multivariate normal distribution as target pdf
pdf <- function(x) dmvnorm(x, c(0,0), matrix(c(1,0.8,0.8,1), ncol=2, byrow = T))

# call sampler
nsamp <- 50000
npars <- 2

set.seed(1312)
res <- am_sd(prior, pdf, t=nsamp, d=npars)

# plot
burnin <- nsamp*0.1
res_burn_pars <- res$chain[-c(1:burnin),]
res_burn_dens <- res$density[-c(1:burnin)]

q68 <- quantile(res_burn_dens, 1-0.68)
r_q68 <- which(round(res_burn_dens, 2) == round(q68, 2))
q90 <- quantile(res_burn_dens, 1-0.9)
r_q90 <- which(round(res_burn_dens, 2) == round(q90, 2))
q95 <- quantile(res_burn_dens, 1-0.95)
r_q95 <- which(round(res_burn_dens, 2) == round(q95, 2))

par(mfrow=c(2,2))

plot(res_burn_pars, xlab="x1", ylab="x2", pch=20, col="red", xlim = c(-4,4), ylim=c(-4,4))
points(res_burn_pars[r_q68,], col="green", pch=20)
points(res_burn_pars[r_q90,], col="black", pch=20)
points(res_burn_pars[r_q95,], col="blue", pch=20)

plot(res$chain[,1], ylab="x1", type="l", col="red", ylim=c(-4,4))
plot(res$chain[,2], ylab="x2", type="l", col="red", ylim=c(-4,4))

savePlot("Sect4_Fig4_am_sd.png", type="png")



### example from Fig. 5 ###
set.seed(50)
par(mfrow=c(1,5))

# uniform prior between -20 and 20
prior <- function(N,d) replicate(N, runif(d, -20, 20))

# multivariate normal distribution as target pdf
pdf <- function(x) mixture(x)

# sample using RWM model
sam_rwm <- rwm(prior = prior, pdf = pdf, t = 50000, d = 1)
hist(sam_rwm$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), main = "RWM", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)

# sample using AM model
sam_am <- am(prior = prior, pdf = pdf, t = 50000, d = 1)
hist(sam_am$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), main = "AM", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)

# sample using AM-SD model
sam_amsd <- am(prior = prior, pdf = pdf, t = 50000, d = 1)
hist(sam_amsd$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), main = "AM-SD", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)

# sample using DE-MC model
sam_de <- de_mc(prior = prior, pdf = pdf, nc=10, t = 5000, d = 1)
hist(sam_de$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 50, main="DE-MC", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)
# # traceplot
# plot(1:dim(sam_de$chain)[1], sam_de$chain[,,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
# for (i in 1:dim(sam_de$chain)[3])
#   points(1:dim(sam_de$chain)[1], sam_de$chain[,,i], pch=i, col=i)

# sample using DREAM
pdf <- function(x) log(mixture(x))
sam_dream <- dream(fun = "pdf", par.info = list(initial = "user", val_ini = prior(10,1)), nc=10, t = 5000, d = 1)
hist(sam_dream$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 50, main="DREAM", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)

savePlot("Sect4_Fig5.png", type="png")





### Case studies from Sect. 5 ###

## 5.1: 1-D mixture distribution
set.seed(1312)
res <- dream_parallel(fun = "mixture", lik = 1, par.info = list(initial = "latin", min = -20, max = 20, prior = "flat", bound = NULL),
                  nc=3, t = 5000, d = 1, adapt = 0.5, burnin = 0, updateInterval = 100,
                  past_sample = T, psnooker = 0.1, mt = 1, ncores = 1)

# Plot
par(mfrow=c(4,2), mar=c(4,4,0,0))
hist(res$chain[,1,], probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 80, main="", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)
# traceplot
plot(1:dim(res$chain)[1], res$chain[,1,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1], res$chain[,1,i], pch=i, col=i)
# prior
plot(1:dim(res$chain)[1], res$chain[,"lp",1], type="n", ylab = "P(x)", xlab = "Sample number of chain")
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1],  res$chain[,"lp",i], pch=i, col=i)
# log-lik
plot(1:dim(res$chain)[1], res$chain[,"ll",1], type="n", ylab = "L(x|Y)", xlab = "Sample number of chain")
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1],  res$chain[,"ll",i], pch=i, col=i)
# post
plot(1:dim(res$chain)[1], res$chain[,"lpost",1], type="n", ylab = "P(x)+L(x|Y)", xlab = "Sample number of chain")
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1],  res$chain[,"lpost",i], pch=i, col=i)
# AR
plot(res$AR[,1], res$AR[,2], type="l", ylab="Acceptance rate", xlab="Number of proposal evaluations")
# CR
plot(res$CR[,1], res$CR[,2], type="n", ylab="Crossover probability", xlab="Number of proposal evaluations", ylim = c(0,1))
for(i in 1:ncol(res$CR))
  lines(res$CR[,1], res$CR[,i], col=i+1)

savePlot("Sect5_1_Fig7_DREAM-zs.png", type="png")

# DEBUG
par(mfrow=c(10,1), mar=c(2,5,0,0))
sub <- 1:5000
# plot(sub, res$chain[sub,,1], type="n", ylab="x", ylim = c(-25,25))
# for (i in 1:dim(res$chain)[3])
#   points(sub, res$chain[sub,,i], pch=i, col=i)
plot(sub, res$DEBUG$gamma[sub,1], type="n", ylab="gamma")
for (i in 1:10)
  points(sub, res$DEBUG$gamma[sub,i], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$lambda[sub,1,1], type="n", ylab="lambda")
for (i in 1:10)
  points(sub, res$DEBUG$lambda[sub,i,1], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$zeta[sub,1,1], type="n", ylab="zeta")
for (i in 1:10)
  points(sub, res$DEBUG$zeta[sub,i,1], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$jump_diff[sub,1,1], type="n", ylab="jump diff")
for (i in 1:10)
  points(sub, res$DEBUG$jump_diff[sub,i,1], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$dx[sub,1,1], type="n", ylab="dx (proposal)")
for (i in 1:10)
  points(sub, res$DEBUG$dx[sub,i,1], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$dx[sub,1,1], type="n", ylab="dx (effective)")
for (i in 1:10)
  points(sub, res$DEBUG$dx_eff[sub,i,1], pch=19, cex=0.5, col=i)
plot(sub, res$DEBUG$std[sub,1], type="l", ylim=c(0, max(res$DEBUG$std[sub,])), ylab="sd of x over chains")
plot(sub, res$DEBUG$J[sub,1], type="n", ylim=c(0, max(res$DEBUG$J[sub,])), ylab="J")
for (i in 1:ncol(res$DEBUG$J))
  lines(sub, res$DEBUG$J[sub,i], col=i)
plot(sub, res$CR[sub,1], type="n", ylab="Crossover probability", ylim = c(0,1))
for(i in 1:ncol(res$CR))
  lines(sub, res$CR[sub,i], col=i+1)
plot(sub, apply(res$AR[sub,], 1, mean), type="l", ylab="Acceptance rate", ylim = c(0,.5))
for(i in 1:ncol(res$AR))
  lines(sub, res$AR[sub,i], col=i+1)

savePlot("Sect5_1_Fig7_detailed.png", type="png")


# compare to BayesianTools package
library(BayesianTools)
prior <- createUniformPrior(-20, 20)
lik <- function(x) log(mixture(x))
bs <- createBayesianSetup(likelihood = lik, prior = prior)
set.seed(1312)
res_BT <- DREAM(bs, settings = list(iterations = 50000,
                                    nCR = 3,
                                    eps = 1e-12, # c_star
                                    pCRupdate = T, updateInterval = 10,
                                    burnin = 0, adaptation = 0.1, # burnin
                                    e = 0.05, # lambda?! In dream() randomly sampled at each iteration
                                    DEpairs = 3, # delta?! In dream() randomly sampled at each iteration
                                    startValue = matrix(replicate(10, runif(1, -20, 20)), ncol=1)
))
res_BT_par <- sapply(res_BT$chains, function(x) x[,"par 1"])
# plot
par(mfrow=c(1,2))
hist(res_BT_par, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 80, main="", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)
# traceplot
plot(1:dim(res_BT_par)[1], res_BT_par[,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res_BT_par)[2])
  points(1:dim(res_BT_par)[1], res_BT_par[,i], pch=i, col=i)

savePlot("Sect5_1_Fig7_BayesianTools.png", type="png")

# DEBUG
par(mfrow=c(5,1), mar=c(2,5,0,0))
plot(1:dim(res_BT_par)[1], res_BT_par[,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res_BT_par)[2])
  points(1:dim(res_BT_par)[1], res_BT_par[,i], pch=i, col=i)
plot(1:dim(res_BT_par)[1], res_BT$gamma[,1], type="n", ylab="gamma")
for (i in 1:10)
  points(1:dim(res_BT_par)[1], res_BT$gamma[,i], pch=19, cex=0.5, col=i)
plot(1:dim(res_BT_par)[1], res_BT$zeta[,1,1], type="n", ylab="zeta")
for (i in 1:10)
  points(1:dim(res_BT_par)[1], res_BT$zeta[,i,1], pch=19, cex=0.5, col=i)
plot(1:dim(res_BT_par)[1], res_BT$jump_diff[,1,1], type="n", ylab="jump diff")
for (i in 1:10)
  points(1:dim(res_BT_par)[1], res_BT$jump_diff[,i,1], pch=19, cex=0.5, col=i)
plot(1:dim(res_BT_par)[1], res_BT$out_cr[,1], type="n", ylab="CR used")
for (i in 1:10)
  points(1:dim(res_BT_par)[1], res_BT$out_cr[,i], pch=19, cex=0.5, col=i)

savePlot("Sect5_1_Fig7_BayesianTools_detailed.png", type="png")




## 5.2: 100-D t-distribution
# TODO: not sure if t_distribution() is implemented correctly, compare with Matlab code!
# set.seed(1312)
# res <- dream(fun = "t_distribution", par.info = list(initial = "latin", min = rep(-5, 100), max = rep(15, 100)),
#              nc=50, t = 100, d = 100, adapt = 0, burnin = 0, updateInterval = 10, thin=5)
#
# # PLOT
# layout(matrix(c(1,2,3,4,5,5), 2, 3, byrow = TRUE), respect = T)
# for (i in c(25,50,75,100)) {
#   # marginal distribution of parameter
#   hist(res$chain[,i,], probability = T, col="red", xlim = c(-4,4), ylim = c(0, 0.4), nclass = 100, xlab=paste0("x", i), main="")
#   # target function
#   lines(seq(-4,4, 0.1), pdf(matrix(seq(-4,4, 0.1), ncol=1)), lwd=2)
# }
# # Gelman-Rubin
# plot(c(1:nrow(res$R_stat))*thin*nc, res$R_stat[,1], type="n", ylab="R-statistic", xlab="No. of function evaluations")
# abline(h=1.2, lty=2, lwd=3)
# for (i in 1:ncol(res$R_stat))
#   lines(c(1:nrow(res$R_stat))*thin*nc, res$R_stat[,i], col=i)





### 25-D trimodal target (Fig. 3, Laloy and Vrugt, 2012) ###
# NOTE: results depend strongly on certain algorithmic parameter selections! Otherwise qualitatively not satisfying results might be produced.
library(mvtnorm)
# func <- function(x) {
#   # 1/3 * dmvnorm(x, rep(-5, length(x)), diag(length(x))) +
#   #   2/3 * dmvnorm(x, rep(5, length(x)), diag(length(x)))
#   #max(log(.Machine$double.xmin), log(
#   3/6 * dmvnorm(x, rep(15, length(x)), diag(length(x))) +
#     2/6 * dmvnorm(x, rep(5, length(x)), diag(length(x))) +
#     1/6 * dmvnorm(x, rep(-5, length(x)), diag(length(x)))
# }
func <- function(x) {
  # NOTE: log-ll needs to be computed directly (instead of just the likelihood) as numerical underflow generates too many zero values which prevents adequate sampling!
  x_t <- t(x)
  mu <- matrix(rep(c(15,5,-5), length(x)), nrow=3)
  lts <- matrix(rep(var2LTsigma(diag(length(x))), 3), nrow=3, byrow=T)
  logL(x_t, pi = c(3/6,2/6,1/6), Mu = mu, LTSigma = lts)
}
set.seed(1312)

res <- dream_parallel(fun = "func", lik = 2,
                      par.info = list(initial = "normal", mu=rep(0,10), cov=5*diag(10), prior = "flat", bound = NULL),
                      nc=5, t = 10000, d = 10, adapt = 0.2, burnin = 0, updateInterval = 100,
                      past_sample = T, archive_update = 10, m0 = 100, psnooker = 0.1, mt=1, ncores = 1)

par(mfrow=c(3,2), mar=c(4,4,0,0))
hist(res$chain[,2,], probability = T, col="red", xlim = c(-10,20), ylim = c(0, 0.3), nclass = 80, main="", xlab="x")
lines(seq(-10, 20, 0.1), exp(sapply(seq(-10, 20, 0.1), func)), lwd=2)
plot(1:dim(res$chain)[1], res$chain[,1,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1], res$chain[,1,i], pch=i, col=i)
# # prior
# plot(1:dim(res$chain)[1], res$chain[,"lp",1], type="n", ylab = "P(x)", xlab = "Sample number of chain")
# for (i in 1:dim(res$chain)[3])
#   points(1:dim(res$chain)[1],  res$chain[,"lp",i], pch=i, col=i)
# log-lik
plot(1:dim(res$chain)[1], res$chain[,"ll",1], type="n", ylab = "L(x|Y)", xlab = "Sample number of chain")
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1],  res$chain[,"ll",i], pch=i, col=i)
# # post
# plot(1:dim(res$chain)[1], res$chain[,"lpost",1], type="n", ylab = "P(x)+L(x|Y)", xlab = "Sample number of chain")
# for (i in 1:dim(res$chain)[3])
#   points(1:dim(res$chain)[1],  res$chain[,"lpost",i], pch=i, col=i)
# AR
plot(res$AR[,1], res$AR[,2], type="l", ylab="Acceptance rate", xlab="Number of proposal evaluations")
# CR
plot(res$CR[,1], res$CR[,2], type="n", ylab="Crossover probability", xlab="Number of proposal evaluations", ylim = c(0,1))
for(i in 1:ncol(res$CR))
  lines(res$CR[,1], res$CR[,i], col=i+1)



# compare to BayesianTools package
library(BayesianTools)
prior <- 0
lik <- function(x) log(func(x))
bs <- createBayesianSetup(likelihood = lik, prior = NULL, parallel=F, lower = rep(-10, 25), upper = rep(20,25))
set.seed(1312)
res_BT <- DREAMzs(bs, settings = list(iterations = 50000,
                                    # nCR = 3,
                                    # eps = 1e-12, # c_star
                                    pCRupdate = T, updateInterval = 10,
                                    burnin = 0, adaptation = 0.2, # burnin
                                    # e = 0.05, # lambda?! In dream() randomly sampled at each iteration
                                    DEpairs = 3, # delta?! In dream() randomly sampled at each iteration
                                    startValue = 5,
                                    ZupdateFrequency = 10,
                                    pSnooker = 0.1,
                                    parallel=F
))
# res_BT <- DREAM(bs, settings = list(iterations = 50000,
#                                     nCR = 3,
#                                     eps = 1e-12, # c_star
#                                     pCRupdate = T, updateInterval = 10,
#                                     burnin = 0, adaptation = 0.2, # burnin
#                                     e = 0.05, # lambda?! In dream() randomly sampled at each iteration
#                                     DEpairs = 3, # delta?! In dream() randomly sampled at each iteration
#                                     startValue = 25
#                                     #startValue = matrix(replicate(10, runif(1, -20, 20)), ncol=1)
# ))
res_BT_par <- sapply(res_BT$chains, function(x) x[,"par 1"])
# plot
par(mfrow=c(1,2))
hist(res_BT_par, probability = T, col="red", xlim = c(-10,20), ylim = c(0, 0.3), nclass = 80, main="", xlab="x")
# plot target function
lines(seq(-10, 20, 0.1), sapply(seq(-10, 20, 0.1), func), lwd=2)
# traceplot
plot(1:dim(res_BT_par)[1], res_BT_par[,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res_BT_par)[2])
  points(1:dim(res_BT_par)[1], res_BT_par[,i], pch=i, col=i)


# dream package
library(dream)
res <- dream::dream(func, func.type = "posterior.density", pars = list(lower = rep(-10, 25), upper = rep(20,25)),
                    INIT = LHSInit, INIT.pars = list(pars = 25, nseq =25),
                    control=list(ndim = 25, nseq=25, DEpairs = 3, nCR=3, boundHandling = "none", burnin.length = 0.2))
res_par <- res_BT_par <- sapply(res_BT$chains, function(x) x[,"par 1"])







### 200-D normal distribution (Fig. 2, Laloy and Vrugt, 2012) ###
library(mvtnorm)
func <- function(x) dmvnorm(x, rep(0, length(x)), diag(1:length(x)) + 0.5 - diag(0.5,length(x),length(x)), log = T)

set.seed(1312)

# DREAM
# res <- dream(fun = "func", lik = 1, par.info = list(initial = "latin", min = -5, max = 15, prior = "uniform", bound = NULL),
#              nc=25, t = 1000, d = 200, adapt = 0.2, burnin = 0, updateInterval = 10)

res <- dream_parallel(fun = "func", lik = 2,
                         par.info = list(initial = "latin", min = -5, max = 15, prior = "uniform", bound = NULL),
                         nc=3, t = 100000, d = 200, adapt = 0.2, burnin = 0, updateInterval = 10,
                         past_sample = T, archive_update = NULL, m0 = NULL, psnooker = 0.1, mt=1, ncores = 1)

par(mfrow=c(1,2))
hist(res$chain[,1,], probability = T, col="red", xlim = c(-10,10), ylim = c(0, 0.5), nclass = 80, main="", xlab="x")
lines(seq(-5, 5, 0.1), sapply(seq(-5, 5, 0.1), function(x) dnorm(x, 0, 1)), lwd=2)
hist(res$chain[,200,], probability = T, col="red", xlim = c(-50,50), ylim = c(0, 0.2), nclass = 80, main="", xlab="x")
lines(seq(-50, 50, 0.1), sapply(seq(-50, 50, 0.1), function(x) dnorm(x, 0, sqrt(200))), lwd=2)

library(BayesianTools)
prior <- createUniformPrior(rep(-5, 200), rep(15,200))
bs <- createBayesianSetup(likelihood = func, prior = prior, parallel=F)
set.seed(1312)
res_BT <- DREAMzs(bs, settings = list(iterations = 50000,
                                      nCR = 3,
                                      eps = 1e-12,
                                      pCRupdate = T, updateInterval = 10,
                                      burnin = 0, adaptation = 0.2,
                                      e = 0.05,
                                      DEpairs = 3,
                                      startValue = 10,
                                      ZupdateFrequency = 10,
                                      pSnooker = 0.1
))
res_BT <- DREAM(bs, settings = list(iterations = 50000,
                                    nCR = 3,
                                    eps = 1e-12, # c_star
                                    pCRupdate = T, updateInterval = 10,
                                    burnin = 0, adaptation = 0.2, # burnin
                                    e = 0.05, # lambda?! In dream() randomly sampled at each iteration
                                    DEpairs = 3, # delta?! In dream() randomly sampled at each iteration
                                    startValue = 25
                                    #startValue = matrix(replicate(10, runif(1, -20, 20)), ncol=1)
))






### Test simple hydrological model ###
fun_wrap <- function(x, forc, obs) {
  # evaluate hydro model
  res <- HydroModel(forcing = forc, param = x["K"])
  if(any(is.na(res)))
    stop("NA in res of fun detected!")
  # calculate log-likelihood of the observations given the model results
  lkh <- sum( dnorm(obs, mean = res, sd = x["sigma"], log = TRUE) )
  # calculate posterior log-pdf ~ lkh (for non-informative prior)
  post <- lkh
  # return in normal space
  return(list(lp = post, sim = res)) # keep_sim = T
  #return(post) # keep_sim = F
}
library(hydrostats)
calc_signatures <- function(dis, prec) {
  # remove (close to) zero flows (for some calculations)
  disnozero <- dis[!(dis <= 0.1)]
  # data.frame object required for external hydrostats-functions; use fake dates
  df.hydrostats <- data.frame(Q=as.numeric(dis), Date=as.Date(1:length(dis), origin=as.Date("1950-01-01")))
  # calculate hydrostats-based statistics (relevant information assigned later on)
  low.hydrostats <- low.spells(df.hydrostats, plot=F)
  high.hydrostats <- high.spells(df.hydrostats, plot=F)
  # MAGNITUDE
  # runoff ratio over the whole given period
  rc <- sum(dis) / sum(prec)
  # exceedance probability of (close to) zero flows
  p_z <- length(disnozero)/length(dis)
  # high flow: average annual maximum flow [V/T]
  # m_h <- high.hydrostats$avg.max.ann
  # FLOW REGIME
  # average slope of flow duration curve for medium range (33% to 66%) of non-zero flows
  # the higher the value the more variable the flow regime; see Sawicz et al. (2011)
  quants <- quantile(dis, probs=1-c(0.33,0.66))
  fdc_s <- (log(quants[1]) - log(quants[2])) / (0.66-0.33)
  # FREQUENCY
  # low flow: average number of low flow events per year [1/year]
  f_l <- low.hydrostats$low.spell.freq
  # high flow: average number of high spell events per year [1/year]
  f_h <- high.hydrostats$spell.freq
    # CONCENTRATION
  # Rate of change in flow events
  # averge absolute daily change during rise periods within high flow events [V/T/day]
  r_r <- high.hydrostats$avg.rise
  r_f <- high.hydrostats$avg.fall
  return(c(rc, p_z, fdc_s, f_l, f_h, r_r, r_f))
}
fun_wrap_abc <- function(x, forc) {
  # evaluate hydro model
  res <- HydroModel(forcing = forc, param = x["K"])
  if(any(is.na(res)))
    stop("NA in res of fun detected!")
    # calculate diagnostic values
  sign <- calc_signatures(res, forc)
  # return
  return(sign)
}
abc_distfun <- function(sim, obs) abs(sim - obs)
fun_wrap_glue <- function(x, forc) {
  # evaluate hydro model
  res <- HydroModel(forcing = forc, param = x["K"])
  if(any(is.na(res)))
    stop("NA in res of fun detected!")
  # return
  return(res)
}
fun_lik <- function(sim, obs) {
  n <- length(sim)
  out <- - n/2 * log( sum((sim-obs)^2) ) # SLS approach with sigma = 1
  return(out)
}

# Read Data
dat=read.table('/home/tobias/Promotion/Workshops_DR/2016/Hydrocourse_Luxemburg/material/Day_2/Exercises/MaimaiCatchment/Maimai.dat',skip=7)
calib=105:304 # calibration period
valid=305:505 # validation period
precip=dat[,7] # precipitation (mm/d)
q_meas=dat[,9] # discharge (mm/d)

d <- 1
nc <- 10
t <- 1000
burnin=0
#obs <- calc_signatures(q_meas[calib], precip[calib])
obs <- q_meas[calib]
set.seed(1312)
res <- dream(fun = "fun_wrap_glue", lik=99, forc = precip[calib],
             par.info = list(initial = "latin", min = 0.0001, max = 1000, mu=1,
                             cov=0.1, bound = "bound", names = "K",
                             prior = "uniform"),
             nc = nc, t = t, d = d, burnin = burnin, adapt = 0.1, beta0 = 1, thin = 1,
             obs = obs, abc_rho = "abc_distfun",
             #abc_e = c(0.3, 0.3, 2, 2, 2, 3, 3) # lik = 22, represents tolerance deviations around obs in which the simulations should lie
             #abc_e = c(0.1, 0.1, 1, 1, 1, 3, 2) # lik = 21, represents sd (measurement error) of streamflow signatures (assumed to be normally distributed)
             #glue_shape = 50 # lik = 31; small values result in high parameter uncertainty; high values produce narrow posterior pdfs of parameters
             lik_fun = "fun_lik" # lik = 99, self-defined likelihood function
             )

# res <- dream(fun = "fun_wrap", forc = precip[calib], obs = q_meas[calib],
#              par.info = list(initial = "normal", min = rep(0.0001,2), max = rep(1000,2), mu=c(1,5), cov=matrix(c(0.1,0,0,0.5), ncol=2, byrow = T), bound = "bound", names = c("K", "sigma")),
#              nc = nc, t = t, d = d, burnin = burnin, adapt = 0.1, beta0 = 1, thin = 1, keep_sim = T)

# traceplots
par(mfrow = c(d, 2))
plot_iter <- 101:nrow(res$chain)
for(p in 1:d) {
  plot(plot_iter, res$chain[plot_iter,p,1], type="n", xlab="Iteration", ylab=paste("par", p), ylim=c(0,20))
  for(i in 1:nc)
    points(plot_iter, res$chain[plot_iter,p,i], pch=19, cex=.5, col=i)
}
# marginal distributions
for(p in 1:d)
  hist(res$chain[plot_iter,p,], probability = T, col="red", nclass = 80, xlab="x", main=paste("par", p))

savePlot("HydroModel_DREAM_trace_margdist.png", type="png")

# histogramms of summary statistics (ABC approach)
par(mfrow = c(1, dim(res$fx)[3]))
for (i in 1:dim(res$fx)[3])
  hist(res$fx[plot_iter,,i], probability = T, col="red", main="")

# acceptability
par(mfrow=c(1,1))
plot(res$AR[,1], type="n")
for(i in 1:3)
  lines(res$AR[,i], col=i)

# get MAP
map_ind <- which(res$chain[,"lpost",] == max(res$chain[,"lpost",]), arr.ind = T)
par_map <- res$chain[map_ind[1,"row"], 1:d,map_ind[1,"col"]]
# simulations with MAP parameters for validation period
sim_map <- HydroModel(precip[valid], par_map[1])
# Plot simulations: obs, sim with MAP and 90% uncertainty bounds of parametric and total uncertainty
#pars_mat <- cbind(c(res$chain[,1,]), c(res$chain[,2,]))
pars_mat <- cbind(c(res$chain[,1,]), NULL)
paramU <- matrix(NA, nrow=nrow(pars_mat), ncol=length(valid))
#totalU <- matrix(NA, nrow=nrow(pars_mat), ncol=length(valid))
for (i in 1:nrow(pars_mat)) {
  # plot the simulation corresponding to the ith parameter set
  K=pars_mat[i,1]
  #sigma=pars_mat[i,2]
  Ysim=HydroModel(precip[valid],K)
  #noise=rnorm(length(valid), sd=sigma)
  paramU[i,] <- Ysim
  #totalU[i,] <- Ysim + noise
}
#total_lo=apply(totalU, 2, quantile, probs=0.05)
#total_hi=apply(totalU, 2, quantile, probs=0.95)
param_lo=apply(paramU, 2, quantile, probs=0.05)
param_hi=apply(paramU, 2, quantile, probs=0.95)
par(mfrow=c(1,1))
plot(valid, sim_map, type="n", ylab="Target variable", xlab="Time step", ylim=c(0, max(q_meas[valid], param_hi)))
#polygon(c(valid, rev(valid)), c(total_lo, rev(total_hi)), col="yellow", border="yellow")
polygon(c(valid, rev(valid)), c(param_lo, rev(param_hi)), col="green", border="green")
lines(valid, sim_map, col="black")
lines(valid, q_meas[valid], col="red")
legend("topleft", legend=c("obs", "sim mode", "param unc", "total unc"), col=c("red", "black", "green", "yellow"), lty=1)

savePlot("HydroModel_DREAM_predict_uncertainty.png", type="png")


# BayesianTools
prior_dens <- function(x) {
  p1 <- dunif(x[1], min = 0.0001, max = 1000, log = T)
  p2 <- dunif(x[2], min = 0.0001, max = 1000, log = T)
  return(p1 + p2)
}
prior_samp <- function(n=1) {
  p1 <- runif(n, min = 0.0001, max = 1000)
  p2 <- runif(n, min = 0.0001, max = 1000)
  return(cbind(p1,p2))
}
lik <- function(x) {
  # evaluate hydro model
  res <- HydroModel(forcing = precip[calib], param = x[1])
  # calculate log-likelihood of the observations given the model results
  lkh <- sum( dnorm(q_meas[calib], mean = res, sd = x[2], log = TRUE) )
  # return
  return(lkh)
}
prior <- createPrior(density = prior_dens, sampler = prior_samp, lower = 0.0001, upper = 1000)
bs <- createBayesianSetup(likelihood = lik, prior = prior)
#set.seed(1312)
res_BT <- DREAM(bs, settings = list(iterations = t*nc,
                                    f = 2.38,
                                    eps = 1e-12, # c_star
                                    burnin = 0,
                                    startValue = mvrnorm(nc, mu=c(1,5), Sigma = matrix(c(0.1,0,0,0.5), ncol=2, byrow = T))
))

# traceplots
par(mfrow = c(d, 1))
res_BT_par_all <- NULL
plot_iter <- 1:t
for(p in 1:d) {
  res_BT_par <- sapply(res_BT$chains, function(x) x[plot_iter,paste("par", p)])
  res_BT_par_all_t <- c(res_BT_par)
  res_BT_par_all <- cbind(res_BT_par_all, res_BT_par_all_t)
  plot(plot_iter, res_BT_par[,1], type="n", xlab="Iteration", ylab=paste("par", p), ylim=c(0,10))
  for(i in 1:nc)
    points(plot_iter, res_BT_par[,i], pch=19, cex=.5, col=i)
}
# parameter plot
plot(res_BT_par_all[-c(1:100),])
