
# Script reproducing examples from Vrugt (2016) Environ. Modell. Softw.

library(MASS)
library(mvtnorm)
#library(lhs)
library(HydroBayes)

setwd("/home/tobias/R/MyPackages/HydroBayes/examples")


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
pdf <- function(x) log(mixture(x))
res <- dream(fun = "pdf", par.info = list(initial = "latin", min = -20, max = 20),
                   nc=10, t = 5000, d = 1, adapt = 0.1, burnin = 0, updateInterval = 10, DEBUG = TRUE)

# Plot
par(mfrow=c(1,2))
hist(res$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 80, main="", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)
# traceplot
plot(1:dim(res$chain)[1], res$chain[,,1], type="n", ylab="x", xlab="Sample of Markov chain", ylim = c(-25,25))
for (i in 1:dim(res$chain)[3])
  points(1:dim(res$chain)[1], res$chain[,,i], pch=i, col=i)

savePlot("Sect5_1_Fig7.png", type="png")

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





## Test simple hydrological model
fun_wrap <- function(x, forc, obs) {
  # evaluate hydro model
  res <- HydroModel(forcing = forc, param = x[1])
  if(any(is.na(res)))
    stop("NA in res of fun detected!")
  # calculate log-likelihood of the observations given the model results
  lkh <- sum( dnorm(obs, mean = res, sd = x[2], log = TRUE) )
  # calculate posterior log-pdf ~ lkh (for non-informative prior)
  post <- lkh
  # return in normal space
  return(list(lp = post, sim = res)) # keep_sim = T
  #return(post) # keep_sim = F
}

# Read Data
dat=read.table('/home/tobias/Promotion/Workshops_DR/2016/Hydrocourse_Luxemburg/material/Day_2/Exercises/MaimaiCatchment/Maimai.dat',skip=7)
calib=105:304 # calibration period
valid=305:505 # validation period
precip=dat[,7] # precipitation
q_meas=dat[,9] # discharge

d <- 2
nc <- 10
t <- 1000
burnin=0.1
set.seed(1312)
res <- dream(fun = "fun_wrap", forc = precip[calib], obs = q_meas[calib],
             par.info = list(initial = "normal", min = rep(0.0001,2), max = rep(1000,2), mu=c(1,5), cov=matrix(c(0.1,0,0,0.5), ncol=2, byrow = T), bound = "bound"),
             nc = nc, t = t, d = d, burnin = burnin, adapt = 0.1, beta0 = 1, thin = 1, keep_sim = T)

# traceplots
par(mfrow = c(d, 2))
plot_iter <- 1:nrow(res$chain)
for(p in 1:d) {
  plot(plot_iter, res$chain[plot_iter,p,1], type="n", xlab="Iteration", ylab=paste("par", p), ylim=c(0,20))
  for(i in 1:nc)
    points(plot_iter, res$chain[plot_iter,p,i], pch=19, cex=.5, col=i)
}
# marginal distributions
for(p in 1:d)
  hist(res$chain[plot_iter,p,], probability = F, col="red", nclass = 80, xlab="x", main=paste("par", p))

savePlot("HydroModel_DREAM_trace_margdist.png", type="png")


# get MAP
map_ind <- which(res$density == max(res$density), arr.ind = T)
par_map <- res$chain[map_ind[1,"row"], ,map_ind[1,"col"]]
# simulations with MAP parameters for validation period
sim_map <- HydroModel(precip[valid], par_map[1])
# Plot simulations: obs, sim with MAP and 90% uncertainty bounds of parametric and total uncertainty
pars_mat <- cbind(c(res$chain[,1,]), c(res$chain[,2,]))
paramU <- matrix(NA, nrow=nrow(pars_mat), ncol=length(valid))
totalU <- matrix(NA, nrow=nrow(pars_mat), ncol=length(valid))
for (i in 1:nrow(pars_mat)) {
  # plot the simulation corresponding to the ith parameter set
  K=pars_mat[i,1]
  sigma=pars_mat[i,2]
  Ysim=HydroModel(precip[valid],K)
  noise=rnorm(length(valid), sd=sigma)
  paramU[i,] <- Ysim
  totalU[i,] <- Ysim + noise
}
total_lo=apply(totalU, 2, quantile, probs=0.05)
total_hi=apply(totalU, 2, quantile, probs=0.95)
param_lo=apply(paramU, 2, quantile, probs=0.05)
param_hi=apply(paramU, 2, quantile, probs=0.95)
par(mfrow=c(1,1))
plot(valid, sim_map, type="n", ylab="Target variable", xlab="Time step")
polygon(c(valid, rev(valid)), c(total_lo, rev(total_hi)), col="yellow", border="yellow")
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
