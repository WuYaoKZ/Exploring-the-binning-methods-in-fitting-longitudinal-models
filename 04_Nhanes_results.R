############ Subsample NHANES multi-level data and fit four types of binning methods ############
# Use a fixed bin size (B=10) when fitting binning methods to NHANES data

rm(list=ls())
library(tidyverse)
library(refund)
library(dplyr)
library(fastFMM)
library(lme4)
set.seed(1234)



# Read in NHANES multi-level data ----

new_nhanes <- readRDS("./data/nhanes_fda_with_r_ml.rds")



# Focus on the young adults ----
# (age between 18 and 30 years old)

select_idx <- new_nhanes$age>=18 & new_nhanes$age<=30
cut_nhanes <- new_nhanes[select_idx,]
cut_nhanes$dayofweek <- factor(cut_nhanes$dayofweek, 
                               levels = c("1","2","3","4","5","6","7"), 
                               labels = c("MON","TUE","WED","THU","FRI","SAT","SUN"))



## Fit four methods (3 binning methods + 1 unbin method) on the sub-sampled data ----

# Covariates: day of wear, gender, age, day of week


### Unbin method ----

ptm <- proc.time()
fit1 <- fui(MIMS ~ dayofwear + gender + age + dayofweek + (1|SEQN), data = cut_nhanes, 
            family = "gaussian", analytic = FALSE, num_boots = 50)  # use bootstrap method to get the confidence band
time1 <- (proc.time() - ptm)[3]

plot1 <- plot_fui(fit1, num_row = 2, return = TRUE)



### Bin average method ----

L <- ncol(cut_nhanes$MIMS)  # 1440 minutes
bin_size <- 10
L_bin <- L/bin_size  # 144

bin_ave_nhanes <- cut_nhanes
bin_ave_nhanes$bin_ave <- matrix(NA, nrow(bin_ave_nhanes), L_bin)  # binned and averaged MIMS

ptm <- proc.time()
for(j in 1:L_bin){
  bin_ave_nhanes$bin_ave[,j] <- rowMeans(bin_ave_nhanes$MIMS2[,(bin_size*(j-1)+1):(bin_size*j)])
}

print("Model: bin and average")    
fit2 <- fui(bin_ave ~ dayofwear + gender + age + dayofweek + (1|SEQN), data = bin_ave_nhanes, family = "gaussian")
time2 <- (proc.time() - ptm)[3]

plot2 <- plot_fui(fit2, num_row = 2, return = TRUE)



### Bin center method ----

bin_center_nhanes <- cut_nhanes
bin_center_nhanes$bin_center <- matrix(NA, nrow(bin_center_nhanes), L_bin)  # binned and centered MIMS

ptm <- proc.time()
for(i in 1:nrow(bin_center_nhanes)){
  bin_center_nhanes$bin_center[i,] <- bin_center_nhanes$MIMS[i, seq(1:L_bin)*bin_size - floor(bin_size/2)]
}

print("Model: bin center")  
fit3 <- fui(bin_center ~ dayofwear + gender + age + dayofweek + (1|SEQN), data = bin_center_nhanes, family = "gaussian")
time3 <- (proc.time() - ptm)[3]

plot3 <- plot_fui(fit3, num_row = 2, return = TRUE)



### Bin and stack method ----

New_SEQN <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
New_dayofwear <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
New_dayofweek <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
New_gender <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
New_age <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)

New_SEQN <- as.factor(rep(cut_nhanes$SEQN, each = bin_size))  # define covariates for binned and stacked data
New_dayofwear <- rep(cut_nhanes$dayofwear, each = bin_size)
New_dayofweek <- rep(cut_nhanes$dayofweek, each = bin_size)
New_gender <- rep(cut_nhanes$gender, each = bin_size)
New_age <- rep(cut_nhanes$age, each = bin_size)

New_MIMS <- matrix(NA, nrow(cut_nhanes)*bin_size, L_bin)  # binned and stacked MIMS

ptm <- proc.time()
for(i in 1:nrow(cut_nhanes)){
  for(j in 1:L_bin){
    New_MIMS[(bin_size*(i-1)+1):(bin_size*i), j] <- cut_nhanes$MIMS[i, (bin_size*(j-1)+1):(bin_size*j)]
  }
}

bin_stack_nhanes <- data.frame(SEQN=New_SEQN, dayofwear=New_dayofwear, dayofweek=New_dayofweek, gender=New_gender, age=New_age, bin_stack=I(New_MIMS))

print("Model: bin and stack")
fit4 <- fui(bin_stack ~ dayofwear + gender + age + dayofweek + (1|SEQN), data = bin_stack_nhanes, family = "gaussian")
time4 <- (proc.time() - ptm)[3]

plot4 <- plot_fui(fit4, num_row = 2, return = TRUE)



# save the results ----

res_list <- list(fit = list(fit1,fit2,fit3,fit4),       # results from the model
                 time = list(time1,time2,time3,time4),  # model fitting time
                 plot = list(plot1,plot2,plot3,plot4))  # plots for covariate estimates

saveRDS(res_list, "./nhanes_res_data/res_list.rds")


