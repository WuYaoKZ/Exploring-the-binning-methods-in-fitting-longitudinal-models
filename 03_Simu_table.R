########## Generate a summary table for comparing simulation results ############
# Output table results only, actual table created manually

rm(list = ls())
library(dplyr)
library(ggplot2)
set.seed(1234)

# Load simulation results ----

### 1st Row: Changing I ----
load("./simu_res_data/list1_I.RData")
load("./simu_res_data/list2_I.RData")
load("./simu_res_data/list3_I.RData")
load("./simu_res_data/list4_I.RData")
list_I <- list(list1, list2, list3, list4) 

### 2nd Row: Changing J ----
load("./simu_res_data/list1_J.RData")
load("./simu_res_data/list2_J.RData")
load("./simu_res_data/list3_J.RData")
load("./simu_res_data/list4_J.RData")
list_J <- list(list1, list2, list3, list4) 

### 3rd Row: Changing L ----
load("./simu_res_data/list1_L.RData")
load("./simu_res_data/list2_L.RData")
load("./simu_res_data/list3_L.RData")
load("./simu_res_data/list4_L.RData")
list_L <- list(list1, list2, list3, list4) 

### 4th Row: Changing B ----
load("./simu_res_data/list1_B.RData")
load("./simu_res_data/list2_B.RData")
load("./simu_res_data/list3_B.RData")
load("./simu_res_data/list4_B.RData")
list_B <- list(list1, list2, list3, list4) 

# Combine all results ----
list_all <- list(list_I, list_J, list_L, list_B)

results_list <- vector("list", 4)

for(para in 1:length(list_all)){  # loop over 3 testing parameters
  
  list_para <- list_all[[para]]  # results list for specific parameter
  
  results_list[[para]] <- vector("list", 4)
  
  for(set in 1:length(list_para)){  # loop over 4 parameter settings
    
    list_set <- list_para[[set]]
    
    ISE <- list_set[[1]]  # ISE
    ISE_1 <- ISE[[1]]  # unbin
    ISE_2 <- ISE[[2]]  # bin and average
    ISE_3 <- ISE[[3]]  # bin center
    ISE_4 <- ISE[[4]]  # bin and stack
    
    Time <- list_set[[2]]  # Computing Time
    Time_1 <- Time[[1]]
    Time_2 <- Time[[2]]
    Time_3 <- Time[[3]]
    Time_4 <- Time[[4]]
    
    Cover_pw <- list_set[[3]]  # Pointwise Coverage
    Cover_pw_1 <- Cover_pw[[1]]
    Cover_pw_2 <- Cover_pw[[2]]
    Cover_pw_3 <- Cover_pw[[3]]
    Cover_pw_4 <- Cover_pw[[4]]
    
    MISE_1 <- mean(rowMeans(bind_rows(ISE_1)))  # Mean_ISE
    MISE_2 <- mean(rowMeans(bind_rows(ISE_2)))
    MISE_3 <- mean(rowMeans(bind_rows(ISE_3)))
    MISE_4 <- mean(rowMeans(bind_rows(ISE_4)))
    
    PW_cover_1 <- mean(colMeans(Cover_pw_1))  # Pointwise_Coverage
    PW_cover_2 <- mean(colMeans(Cover_pw_2))
    PW_cover_3 <- mean(colMeans(Cover_pw_3))
    PW_cover_4 <- mean(colMeans(Cover_pw_4))
    
    Median_Time_1 <- median(Time_1)  # Median_Time
    Median_Time_2 <- median(Time_2)
    Median_Time_3 <- median(Time_3)
    Median_Time_4 <- median(Time_4)
    
    MISE <- round(c(MISE_1, MISE_2, MISE_3, MISE_4) * 1e05, 2)  # Scale the MISE by multiply 1e05
    PW <- round(c(PW_cover_1, PW_cover_2, PW_cover_3, PW_cover_4), 2)
    MT <- round(c(Median_Time_1, Median_Time_2, Median_Time_3, Median_Time_4), 2)
    
    results_list[[para]][[set]] <- list(MISE, PW, MT)  # outputs: MISE, Point-wise coverage, Median computing time
    
  }
  
}

#str(results_list)



# Results ----

# Three rows for each case: MISE, Coverage Probability, Computing Time


### Table 1: Changing I ----

results_list[[1]][[1]]  # I=50
results_list[[1]][[2]]  # I=100
results_list[[1]][[3]]  # I=200
results_list[[1]][[4]]  # I=400

### Table 2: Changing J ----

results_list[[2]][[1]]  # J=5
results_list[[2]][[2]]  # J=10
results_list[[2]][[3]]  # J=20
results_list[[2]][[4]]  # J=40

### Table 3: Changing L ----

results_list[[3]][[1]]  # L=50
results_list[[3]][[2]]  # L=100
results_list[[3]][[3]]  # L=200
results_list[[3]][[4]]  # L=400

### Table 4: Changing B ----

results_list[[4]][[1]]  # B=5
results_list[[4]][[2]]  # B=10
results_list[[4]][[3]]  # B=20
results_list[[4]][[4]]  # B=40





