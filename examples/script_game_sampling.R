
library(HydroBayes)

# d-variate normal distribution, zero mean, dxd covariance matrix
library(mvtnorm)
#func <- function(x) dmvnorm(x, rep(0, length(x)), diag(1:length(x)) + 0.5 - diag(0.5,length(x),length(x)))
func <- function(x) 1/3 * dmvnorm(x, rep(-5, length(x)), diag(length(x))) + 2/3 * dmvnorm(x, rep(5, length(x)), diag(length(x)))

# run dream to sample unnormalised posterior
d <- 20
res <- dream_parallel(fun = "func", lik = 1,
                      par.info = list(initial = "latin", min = -10, max = 10, prior = "flat", bound = NULL),
                      nc=20, t = 8000, d = d, adapt = 0.5, burnin = 0.5, thin = 1, updateInterval = 10, past_sample = T,
                      psnooker = 0, mt=1)

# plot sampling result together with known target distribution
par(mfrow=c(ceiling(4/3), 3))
for(i in c(1,5,10,15, 20)) {
  hist(res$chain[-1,i,], probability = T, col="red", xlim = c(-15,15), ylim = c(0, 0.5), nclass = 80, main="", xlab="x")
  title(main = paste("par", i))
  #lines(seq(-15, 15, 0.1), sapply(seq(-15, 15, 0.1), function(x) dnorm(x, 0, sqrt(i))), lwd=2)
  lines(seq(-15, 15, 0.1), sapply(seq(-15, 15, 0.1), function(x) func(x)), lwd=2)
}

# prepare input data for GAME sampling algorithm
theta <- apply(res$chain[-1,1:d,], 2, c) # matrix of posterior samples
theta_post <- exp(c(res$chain[-1,"lpost",])) # vector of unnormalised posterior densities, q1 in Volpi et al., 2017

library(EMCluster)
# algorithmic parameters
Jmax <- 5
m1 <- 1000
m0 <- 5000
m_sub <- 4000
omega <- 0
ob <- T
ic="bic"

par.info = list(initial = "latin", min = -10, max = 10, prior = "flat", bound = NULL)
lik = 1
fun = "func"
obs = NULL
abc_rho = NULL
abc_e = NULL
glue_shape = NULL
lik_fun = NULL
ncores = 8

replicate(5, {
res_game <- game(theta = theta, theta_post = theta_post, m_sub = m_sub, m0 = m0, m1 = m1, Jmax = Jmax, ic = ic,
                 omega = omega, ob = ob, par.info = par.info, lik = lik, fun = fun, ncores = ncores)
res_game
1/res_game$Z
})





# # BayesianTools variant
# library(BayesianTools)
# prior <- createUniformPrior(rep(-30, 2), rep(30,2))
# lik <- function(x) log(func(x))
# bs <- createBayesianSetup(likelihood = lik, prior = prior, parallel=F)
# res_BT <- DREAM(bs, settings = list(iterations = 40000,
#                                     nCR = 3,
#                                     eps = 1e-12, # c_star
#                                     pCRupdate = T, updateInterval = 10,
#                                     burnin = 0.5, adaptation = 0.5, # burnin
#                                     e = 0.05, # lambda?! In dream() randomly sampled at each iteration
#                                     DEpairs = 3, # delta?! In dream() randomly sampled at each iteration
#                                     startValue = 10
# ))
#
# par(mfrow=c(1,2))
# for(i in 1:2) {
#   res_BT_par <- sapply(res_BT$chains, function(x) x[,paste("par", i)])
#   hist(res_BT_par, probability = T, col="red", xlim = c(-30,30), ylim = c(0, 0.5), nclass = 80, main="", xlab="x")
#   title(main = paste("par", i))
#   lines(seq(-30, 30, 0.1), sapply(seq(-30, 30, 0.1), function(x) dnorm(x, 0, i)), lwd=2)
# }
#
# evidence_BT <- bridgesample(cbind(theta, log(theta_post)), nParams = 10)


# bridgesampling package
library(bridgesampling)
lower <- rep(-15,10)
names(lower) <- colnames(theta)
upper <- rep(15,10)
names(upper) <- colnames(theta)
res_bs <- bridge_sampler(theta, log_posterior = function(x, data) log(func(x)), lb=lower, ub=upper)
