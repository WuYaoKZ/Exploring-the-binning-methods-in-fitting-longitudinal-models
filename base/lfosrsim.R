#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#' 
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param I number of subjects.
#' @param J mean number of observations per subject.
#' @param L number of grid points on the functional domain.
#' @param beta_true true fixed effects functions.
#' @param psi_true true orthonormal functions to generate subject-specific random effects.
#' @param psi2_true true orthonormal functions to generate subject/visit-specific random effects.
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio in Gaussian data generation
#' @export
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

lfosrsim <- function(family = "gaussian", I = 100, J = 10, L = 100, 
                      beta_true, psi_true, psi2_true = NULL, 
                      SNR_B = 1, SNR_sigma = 1){
  library(mvtnorm)
  
  ## generate number of visits for each subject from poisson distribution
  J_subj <- pmax(rpois(I, J), 1)  # Set 0 visits to 1
  
  ## generate fixed effects
  n <- sum(J_subj)  # sum of all visits
  X_des = cbind(1, rnorm(n, 0, 2))  # one row for each visit [1,x]
  fixef <- X_des %*% beta_true  # beta0 + beta1*x
  
  ## generate random effects
  subj <- as.factor(rep(1:I, J_subj))  # replicate the subject ID by # of visits, and factorize them
  Z_des <- model.matrix( ~ 0 + subj)  # create design matrix
  c_true <- rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5)))  ## simulate score function (principal component)
  b_true <- c_true %*% psi_true
  ranef = Z_des %*% b_true  # based on # of individuals (columns)
  if(!is.null(psi2_true)){  ## by default do not add subject-visit random deviation
    c2_true <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))  # based on # of visits (rows)
    ranef <- ranef + c2_true %*% psi2_true
  }
  
  ## generate linear predictors
  ranef <- sd(fixef)/sd(ranef)/SNR_B*ranef  ## adjust for relative importance of random effects
  eta_true <- fixef + ranef
  
  ## generate longitudinal functional data
  Y_obs <- matrix(NA, n, L)   # matrix made of NA's, n (visit) rows, L (grid) columns
  p_true <- plogis(eta_true)  # prob of event: exp(eta_true)/(1+exp(eta_true))
  lam_true <- exp(eta_true)  # lambda for poisson dist
  sd_signal <- sd(eta_true)
  for(i in 1:n){  # visit level
    for(j in 1:L){  # grid level
      if(family == "gaussian"){
        Y_obs[i, j] <- rnorm(1, mean = eta_true[i, j], sd = sd_signal/SNR_sigma)
      }else if(family == "binomial"){
        Y_obs[i, j] <- rbinom(1, 1, p_true[i, j])
      }else if(family == "poisson"){
        Y_obs[i, j] <- rpois(1, lam_true[i, j])
      }
    }
  }
  
  ## combine simulated data
  visit <- rep(1:I, J_subj)  # replicate the subject ID by # of visits
  for(i in 1:I){
    visit[which(subj == i)] <- 1:J_subj[i]  # change the replicates of subject ID to 1 to subject ID
  }
  dat.sim <- data.frame(ID = subj, visit = visit, X = X_des[,2], Y = I(Y_obs), eta = I(eta_true))  
  # I() means to treat object as it is, inhibit conversion
  
  return(dat.sim)
}
