
# Script reproducing examples from Vrugt (2016) Environ. Modell. Softw.

library(MASS)
library(mvtnorm)
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
burnin <- 0
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

savePlot("Fig4_am_sd.png", type="png")



### example from Fig. 5 ###
set.seed(50)
par(mfrow=c(1,4))

# uniform prior between -20 and 20
prior <- function(N,d) t(replicate(N, runif(d, -20, 20)))

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

savePlot("Fig5_test_mixture_distribution.png", type="png")


sam_dream <- dream(prior = prior, pdf = pdf, nc=10, t=5000, d=1)
hist(sam_dream$chain, probability = T, col="red", xlim = c(-12,15), ylim = c(0, 0.4), nclass = 50, main="DREAM", xlab="x")
# plot target function
lines(seq(-12, 15, 0.1), mixture(seq(-12, 15, 0.1)), lwd=2)
