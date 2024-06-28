############ The unified program for getting simulation results for 4 different methods (unbin + binning) ############

# NOTE: This version's program has accounted for the non-integer bin center and inappropriate bin size

# Loading needed packages and r function "lfosrsim" from base folder to simulate longitudinal functional data ----

rm(list = ls())
library(refund)
library(dplyr)
library(fastFMM)
source("./base/lfosrsim.R")
set.seed(1234)

# Write a function to fit four different methods on the simulated data and return the model fitting performance measures ----
# (ISE, computing time, and empirical coverage)

run_simulation <- function(I,J,L,SNR_B,SNR_sigma,family,nsim,bin_size){
  
  ## Check if input parameters are reasonable ----
  
  # (1) check if the bin center is an integer
  if ((bin_size %% 2) == 0) {
    print("Impropriate bin size! Bin center is not an integer!")
    nonint_bc = 1  # flag variable for non-integer bin center
  } else {
    nonint_bc = 0
  }
  cat("check point indicator:", nonint_bc)
  
  # (2) check if the number of bins is integer
  if ((L/bin_size) != round(L/bin_size)) {
    stop("Impropriate bin size! Number of bins is not integer!")
  }
  
  ## Specify true fixed and random effects functions ----
  
  # (1) true fixed effects
  grid <- seq(0, 1, length = L)
  beta_true <- matrix(NA, 2, L)
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
  
  rownames(beta_true) <- c("Intercept", "x")
  
  # (2) true random effects
  psi_true <- matrix(NA, 2, L)
  psi_true[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
  psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
  psi_true[2,] <- sin(4*grid*pi)
  psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))
  
  ## Do simulations on a local laptop ----
  
  # use four lists to store model fitting results on simulated data
  
  # ISE
  sim_ise_unbin <- list()
  sim_ise_bin_average <- list()
  sim_ise_bin_center <- list()
  sim_ise_bin_stack <- list()
  
  # Model fitting time
  sim_time_unbin <- array(NA, length(nsim))
  sim_time_bin_average <- array(NA, length(nsim))
  sim_time_bin_center <- array(NA, length(nsim))
  sim_time_bin_stack <- array(NA, length(nsim))
  
  # Coverage probability
  sim_cover_unbin <- matrix(NA, nrow = nsim, ncol = L) # store pointwise coverage
  sim_cover_bin_average <- matrix(NA, nrow = nsim, ncol = L/bin_size)
  sim_cover_bin_center <- matrix(NA, nrow = nsim, ncol = L/bin_size)
  sim_cover_bin_stack <- matrix(NA, nrow = nsim, ncol = L/bin_size)
  
  for(iter in 1:nsim){
    
    set.seed(iter)
    
    ### Simulate the data for FUI using "lfosrsim" function ----
    data <- lfosrsim(family, I, J, L, beta_true, psi_true, SNR_B = SNR_B, SNR_sigma = SNR_sigma)
    data0 <- data # backu the simulated data
    
    L_bin <- L/bin_size # Dimension of binned data "L_bin"
    
    #### (1) Bin and stack Method ----
    
    ptm <- proc.time()
    
    NewID <- matrix(NA, dim(data)[1]*bin_size, 1) # ID for stacked data
    NewID <- as.factor(rep(data$ID, each = bin_size))
    New_visit <- matrix(NA, dim(data)[1]*bin_size, 1) # Number of visits for stacked data
    New_visit <- rep(data$visit, each = bin_size)
    NewX <- matrix(NA, dim(data)[1]*bin_size, 1) # Fixed effects for stacked data
    NewX <- rep(data$X, each = bin_size)
    NewY <- matrix(NA, dim(data)[1]*bin_size, L_bin) # Functional response for stacked data
    New_eta <- matrix(NA, dim(data)[1]*bin_size, L_bin) # Linear predictor for stacked data
    
    for(i in 1:dim(data)[1]){
      for(j in 1:L_bin){
        NewY[(bin_size*(i-1)+1):(bin_size*i), j] <- data$Y[i, (bin_size*(j-1)+1):(bin_size*j)]
        New_eta[(bin_size*(i-1)+1):(bin_size*i), j] <- data$eta[i, (bin_size*(j-1)+1):(bin_size*j)]
      }
    }
    
    # Create new binned and stacked data
    new_data <- data.frame(NewID=NewID, New_visit=New_visit, NewX=NewX, NewY=I(NewY), New_eta=I(New_eta))
    
    print("Model: bin and stack")
    fit <- fui(NewY ~ NewX + (1 | NewID), data = new_data) # apply FUI approach on binned and stacked data
    time <- (proc.time() - ptm)[3] # computing time
    
    bin_beta <- matrix(NA, dim(beta_true)[1], L_bin)
    
    for(i in 1:dim(beta_true)[1]){
      if(nonint_bc == 0){
        bin_beta[i,] <- beta_true[i, seq(1:L_bin)*bin_size - floor(bin_size/2)] # center point of the bin
      } else  
      if(nonint_bc == 1){
        bin_beta[i,] <- apply(cbind(beta_true[i, seq(1:L_bin)*bin_size - bin_size/2],
                                    beta_true[i, seq(1:L_bin)*bin_size - bin_size/2 + 1]),
                              1, mean) # mean of the two center points of the bin
      }
    }
    
    ISE <- (fit$betaHat - bin_beta)^2  # ISE: integrated squared error
    
    cover_pw <- matrix(FALSE, nrow(bin_beta), L_bin) # coverage probability of pointwise confidence bands
    
    for(p in 1:nrow(bin_beta)){
      cover_upper_pw <- which(fit$betaHat[p,]+1.96*sqrt(diag(fit$betaHat.var[,,p])) > bin_beta[p,]) # upper bound coverage
      cover_lower_pw <- which(fit$betaHat[p,]-1.96*sqrt(diag(fit$betaHat.var[,,p])) < bin_beta[p,]) # lower bound coverage
      cover_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE  # check if CI cover the true value: point-wise coverage
    }
    
    sim_ise_bin_stack[[iter]] <- ISE[2,]
    sim_time_bin_stack[iter] <- time
    sim_cover_bin_stack[iter,] <- cover_pw[2,]
    
    #### (2) Bin and average Method ----
    
    ptm <- proc.time()
    
    NewY <- matrix(NA, dim(data)[1], L_bin) # binned and averaged functional response
    
    for(i in 1:dim(data)[1]){
      for(j in 1:L_bin){
        NewY[i,j] <- mean(data$Y[i,(bin_size*(j-1)+1):(bin_size*j)])
      }
    }
    
    data$NewY <- NewY
    
    print("Model: bin and average")    
    fit <- fui(NewY ~ X + (1 | ID), data = data) # apply FUI approach using the binned and averaged response
    time <- (proc.time() - ptm)[3] # computing time
    
    ISE <- (fit$betaHat - bin_beta)^2  # ISE: integrated squared error 
    
    cover_pw <- matrix(FALSE, nrow(bin_beta), L_bin)
    
    for(p in 1:nrow(bin_beta)){
      cover_upper_pw <- which(fit$betaHat[p,]+1.96*sqrt(diag(fit$betaHat.var[,,p])) > bin_beta[p,]) # upper bound coverage
      cover_lower_pw <- which(fit$betaHat[p,]-1.96*sqrt(diag(fit$betaHat.var[,,p])) < bin_beta[p,]) # lower bound coverage
      cover_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE  # check if CI cover the true value: point-wise
    }
    
    sim_ise_bin_average[[iter]] <- ISE[2,]
    sim_time_bin_average[iter] <- time
    sim_cover_bin_average[iter,] <- cover_pw[2,]
    
    
    #### (3) Bin center Method ----
    
    # Take the bin center
    
    ptm <- proc.time()
    
    Y_bin_center <- matrix(NA, dim(data)[1], L_bin) # Binned and centered functional response
    
    for(i in 1:dim(data)[1]){
      if(nonint_bc == 0){
        Y_bin_center[i,] <- data$Y[i, seq(1:L_bin)*bin_size - floor(bin_size/2)] # center point of the bin
      } else
      if(nonint_bc == 1){
        Y_bin_center[i,] <- apply(cbind(data$Y[i, seq(1:L_bin)*bin_size - bin_size/2],
                                        data$Y[i, seq(1:L_bin)*bin_size - bin_size/2 + 1]),
                                  1, mean) # mean of the two center points of the bin
      }
    }
    data$Y_bin_center <- Y_bin_center
    
    print("Model: bin center")  
    fit <- fui(Y_bin_center ~ X + (1 | ID),data = data)
    time <- (proc.time() - ptm)[3]
    
    ISE <- (fit$betaHat - bin_beta)^2  # ISE: integrated squared error
    
    cover_pw <- matrix(FALSE, nrow(bin_beta), L_bin)
    
    for(p in 1:nrow(bin_beta)){
      cover_upper_pw <- which(fit$betaHat[p,]+1.96*sqrt(diag(fit$betaHat.var[,,p])) > bin_beta[p,]) # upper bound coverage
      cover_lower_pw <- which(fit$betaHat[p,]-1.96*sqrt(diag(fit$betaHat.var[,,p])) < bin_beta[p,]) # lower bound coverage
      cover_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE  # check if CI cover the true value: point-wise
    }
    
    sim_ise_bin_center[[iter]] <- ISE[2,]
    sim_time_bin_center[iter] <- time
    sim_cover_bin_center[iter,] <- cover_pw[2,]
    
    #### (4) Unbin Method ----
    
    print("Model: unbin")
    ptm <- proc.time()
    fit <- fui(Y ~ X + (1 | ID),data = data0)
    time <- (proc.time() - ptm)[3]
    
    ISE <- (fit$betaHat - beta_true)^2  # ISE: integrated squared error 
    
    cover_pw <- matrix(FALSE, nrow(beta_true), L)
    
    for(p in 1:nrow(beta_true)){
      cover_upper_pw <- which(fit$betaHat[p,]+1.96*sqrt(diag(fit$betaHat.var[,,p])) > beta_true[p,]) # upper bound coverage
      cover_lower_pw <- which(fit$betaHat[p,]-1.96*sqrt(diag(fit$betaHat.var[,,p])) < beta_true[p,]) # lower bound coverage
      cover_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE  # check if CI cover the true value: point-wise
    }
    
    sim_ise_unbin[[iter]] <- ISE[2,]
    sim_time_unbin[iter] <- time
    sim_cover_unbin[iter,] <- cover_pw[2,]
    
    print(iter)
    
  }
  
  ### Organize model fitting results into lists ----
  
  ISE_list <- list(sim_ise_unbin,
                   sim_ise_bin_average,
                   sim_ise_bin_center,
                   sim_ise_bin_stack)
  
  Time_list <- list(sim_time_unbin,
                    sim_time_bin_average,
                    sim_time_bin_center,
                    sim_time_bin_stack)
  
  Coverpw_list <- list(sim_cover_unbin,
                       sim_cover_bin_average,
                       sim_cover_bin_center,
                       sim_cover_bin_stack)
  
  Total_list <- list(ISE_list, Time_list, Coverpw_list)
  
  return(Total_list)
}

# Run the function for getting model fitting results: ISE, computing time, and coverage probability ----

## Get simulated results for changing the number of subjects (I) ----
# I = 50, 100, 200, 400

list1 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list2 <- run_simulation(I = c(100),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list3 <- run_simulation(I = c(200),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list4 <- run_simulation(I = c(400),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

### Save the simulated data and model results to local files ###
save(list1, file = "./simu_res_data/list1_I.Rdata")
save(list2, file = "./simu_res_data/list2_I.Rdata")
save(list3, file = "./simu_res_data/list3_I.Rdata")
save(list4, file = "./simu_res_data/list4_I.Rdata")



## Get simulated results for changing the mean number of visits per subject (J) ----
# J = 5, 10, 20, 40

list1 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list2 <- run_simulation(I = c(50),
                        J = c(10),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list3 <- run_simulation(I = c(50),
                        J = c(20),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list4 <- run_simulation(I = c(50),
                        J = c(40),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

### Save the simulated data and model results to local files ###
save(list1, file = "./simu_res_data/list1_J.Rdata")
save(list2, file = "./simu_res_data/list2_J.Rdata")
save(list3, file = "./simu_res_data/list3_J.Rdata")
save(list4, file = "./simu_res_data/list4_J.Rdata")



## Get simulated results for changing the dimension on the functional domain (L) ----
# L = 50, 100, 200, 400

list1 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(50),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list2 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(100),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list3 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(200),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list4 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(400),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

### Save the simulated data and model results to local files ###
save(list1, file = "./simu_res_data/list1_L.Rdata")
save(list2, file = "./simu_res_data/list2_L.Rdata")
save(list3, file = "./simu_res_data/list3_L.Rdata")
save(list4, file = "./simu_res_data/list4_L.Rdata")



## Get simulated results for changing Binsize (B) ----

list1 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(400),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(5))

list2 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(400),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(10))

list3 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(400),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(20))

list4 <- run_simulation(I = c(50),
                        J = c(5),
                        L = c(400),
                        SNR_B = c(0.5),
                        SNR_sigma = c(1),
                        family = c("gaussian"),
                        nsim = c(200),
                        bin_size = c(40))

### Save the simulated data and model results to local files ###
save(list1, file = "./simu_res_data/list1_B.Rdata")
save(list2, file = "./simu_res_data/list2_B.Rdata")
save(list3, file = "./simu_res_data/list3_B.Rdata")
save(list4, file = "./simu_res_data/list4_B.Rdata")




