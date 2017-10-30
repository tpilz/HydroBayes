# calibration of the mixture distribution employing the EM algorithm (external R package EMCluster)
em_pmix <- function(theta, theta_post, J, ic) {
  # apply Expectation-Maximization method to calibrate properties of the Gaussian mixture distribution
  em_res <- init.EM(theta, J) # NOTE: LTSigma counts row-wise staring at upper left of the LT matrix

  # calculate information criterium to evaluate the calibration
  if(grepl("bic", ic, ignore.case = T)) {
    # calculate BIC
    ic_val <- em.bic(theta, em_res)

  } else if(grepl("var", ic, ignore.case = T)) {
    # densities of the mixture model for all observations
    pmix_dens <- dmixmvn(theta, em_res)
    # calculate variance of the likelihood ratios
    ic_val <- var(theta_post / pmix_dens)
  } else stop("Wrong choice of Information Criterium to evaluate the EM calibration of the Gaussian mixture distribution. Must be one oc {'BIC' or 'var'}.")

  # output list
  return(list(pmix = em_res,
              ic_val = ic_val,
              J = J))
}
