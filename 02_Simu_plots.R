########## Generate merged plots for simulation results ############

rm(list = ls())
library(dplyr)
library(ggplot2)
set.seed(1234)

source("./00_get_df.R")

# Read in simulation results for changing I, J, L, B ----

load("./simu_res_data/list1_I.RData")
load("./simu_res_data/list2_I.RData")
load("./simu_res_data/list3_I.RData")
load("./simu_res_data/list4_I.RData")
list_I_1 <- list1
list_I_2 <- list2
list_I_3 <- list3
list_I_4 <- list4

load("./simu_res_data/list1_J.RData")
load("./simu_res_data/list2_J.RData")
load("./simu_res_data/list3_J.RData")
load("./simu_res_data/list4_J.RData")
list_J_1 <- list1
list_J_2 <- list2
list_J_3 <- list3
list_J_4 <- list4

load("./simu_res_data/list1_L.RData")
load("./simu_res_data/list2_L.RData")
load("./simu_res_data/list3_L.RData")
load("./simu_res_data/list4_L.RData")
list_L_1 <- list1
list_L_2 <- list2
list_L_3 <- list3
list_L_4 <- list4

load("./simu_res_data/list1_B.RData")
load("./simu_res_data/list2_B.RData")
load("./simu_res_data/list3_B.RData")
load("./simu_res_data/list4_B.RData")
list_B_1 <- list1
list_B_2 <- list2
list_B_3 <- list3
list_B_4 <- list4

list_all_params <- list(
  list(list_I_1, list_I_2, list_I_3, list_I_4),
  list(list_J_1, list_J_2, list_J_3, list_J_4),
  list(list_L_1, list_L_2, list_L_3, list_L_4),
  list(list_B_1, list_B_2, list_B_3, list_B_4)
)

# Structure of "list_all_params":

### (1) 1st layer: four sample size parameters (I, J, L, B)
### (2) 2nd layer: four settings for each sample size parameter (I/L: 50-100-200-400; J/B: 5-10-20-40)
### (3) 3rd layer: three model performance metrics (ISE, computing time, point-wise coverage)
### (4) 4th layer: four model fitting methods ("unbin", "bin and average", "bin center", "bin and stack")
### (5) 5th layer: 200 times of simulations



# Define a function to extract simulation results ----

get_df <- function(list, item){
  
  list <- list[[item]]  # extract specific metrics (ISE, computing time, point-wise coverage)
  
  if(item==1){  # item=1: ISE;
    unbin <- sapply(list[[1]],mean)  # calculate MISE
    bin_average <- sapply(list[[2]],mean)
    bin_center <- sapply(list[[3]],mean)
    bin_stack <- sapply(list[[4]],mean)
  }
  else if(item==2){  # item=2: computing time;
    unbin <- list[[1]]
    bin_average <- list[[2]]
    bin_center <- list[[3]]
    bin_stack <- list[[4]]
  }
  else if(item==3){  # item=3: pointwise coverage;
    unbin <- colMeans(list[[1]])  # calculate coverage probability
    bin_average <- colMeans(list[[2]])
    bin_center <- colMeans(list[[3]])
    bin_stack <- colMeans(list[[4]])
  }
  
  if(item==1|item==2){  
    df <- data.frame(
      Method = factor(
        rep(c("Unbin", "Bin and average", "Bin center", "Bin and stack"), each = length(unbin)),
        levels = c("Unbin", "Bin and average", "Bin center", "Bin and stack")
      ),
      Values = c(unbin, bin_average, bin_center, bin_stack))
  }
  else if(item==3){
    df <- data.frame(
      Method = factor(
        c(rep(c("Unbin"), each = length(unbin)),
          rep(c("Bin and average", "Bin center", "Bin and stack"), each = length(bin_average))
        ),
        levels = c("Unbin", "Bin and average", "Bin center", "Bin and stack")
      ),
      Values = c(unbin, bin_average, bin_center, bin_stack))
  }
  
  return(df)
}



# Organize and plot results for each sample size parameter (I, J, L,B) ----

### Define vectors and lists for 4 parameters ----

params <- c("I", "J", "L", "B")  # sample size parameters
param_labels <- c("I (Number of subjects)",
                  "J (Mean number of observations per subject)",
                  "L (Number of points on the functional domain)",
                  "B (Bin Size) (L = 400)")  # labels for sample size parameters

metrics <- c("ISE", "Time", "Coverage")  # model performance metrics
reps <- list(c(50, 100, 200, 400), c(5,10,20,40), c(50,100,200,400), c(5,10,20,40))  # number of replicates
plots <- list(plot1 = list(), plot2 = list(), plot3 = list(), plot4 = list())  # lists for storing the plots

legend_position <- list("none", "none", "none", "right")  # Merged plot: only show legend for the last parameter
### legend_position <- list("right", "right", "right", "right")  # Separate plots: show legend for each parameter

for (p in 1:length(params)) { # p for 4 parameters (I, J, K, L)
  
  p_list <- list_all_params[[p]]  # extract results for the p-th parameter
  
  ### Extract results using "get_df" function ----
  p_ISE <- lapply(p_list, get_df, 1)
  p_Time <- lapply(p_list, get_df, 2)
  p_Coverage <- lapply(p_list, get_df, 3)
  
  ### Store results in data frames ----
  df_ISE <- data.frame(
    P = rep(reps[[p]], c(nrow(p_ISE[[1]]), nrow(p_ISE[[2]]), nrow(p_ISE[[3]]), nrow(p_ISE[[4]]))),
    rbind(p_ISE[[1]], p_ISE[[2]], p_ISE[[3]], p_ISE[[4]])
  )
  df_Time <- data.frame(
    P = rep(reps[[p]], c(nrow(p_Time[[1]]), nrow(p_Time[[2]]), nrow(p_Time[[3]]), nrow(p_Time[[4]]))),
    rbind(p_Time[[1]], p_Time[[2]], p_Time[[3]], p_Time[[4]])
  )
  df_pwcover <- data.frame(
    P = rep(reps[[p]], c(nrow(p_Coverage[[1]]), nrow(p_Coverage[[2]]), nrow(p_Coverage[[3]]), nrow(p_Coverage[[4]]))),
    rbind(p_Coverage[[1]], p_Coverage[[2]], p_Coverage[[3]], p_Coverage[[4]])
  )
  
  ### Combine ISE and Time in one data frame ----
  combined_df <- cbind(df_ISE, df_Time[,3])
  colnames(combined_df)[3] <- "ISE"
  colnames(combined_df)[4] <- "Time" 
  combined_df$ISE <- combined_df$ISE*1e05  # scale the ISE value
  
  mean_time <- combined_df %>% 
    group_by(P, Method) %>% 
    summarise(mean_time = mean(Time)) %>%  # calculate mean computing time
    ungroup()
  combined_df <- merge(combined_df, mean_time, by=c("P","Method"),all.x = TRUE)

  ### Plot for ISE and computing time ----
  plots[[p]][[1]] = ggplot(combined_df, aes(
      x = P, y = ISE, 
      fill = Method
    )) +
    geom_boxplot(aes(group = interaction(P, Method))) +
    scale_fill_manual(
      values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
    ) +
    scale_x_continuous(breaks = unique(combined_df$P)) +
    # geom_line(aes(x = P, y = 3*log(1+mean_time), color= Method), linewidth=0.7) +
    # scale_color_manual(
    #   values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
    # ) +
    geom_line(aes(x = P, y = 3*log(1+mean_time), color= Method), linewidth=0.7) +
    scale_color_manual(
      values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
      name = "Computing Time"
    ) +
    labs(
      x = param_labels[[p]]
    ) +
    scale_y_continuous(
      # set limits for the Y axis
      limits = c(0,15),
      # Features of the first axis
      name = expression(paste("ISE(", beta[1], "(s))(x", 10^-5, ")")),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*(1/3), name="Log Scaled Computing Time (log (1 + min))")
    ) + 
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # Change x-axis label font size here
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = legend_position[[p]])  # control the position of legends
  
  ### Plot for point-wise coverage ----
  plots[[p]][[2]] = ggplot(df_pwcover, aes(
      x = P, y = Values*100, 
      fill = Method,
      group = interaction(P, Method)
    )) +
    geom_boxplot() +
    scale_fill_manual(
      values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
    ) +
    scale_x_continuous(breaks = unique(df_pwcover$P)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(40,100), breaks = c(0,25,50,75,95,100)) +
    labs(
      x = param_labels[[p]],
      y = "Pointwise coverage (%)"
    ) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # Change x-axis label font size here
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = legend_position[[p]])  # control the position of legends

}

plots


# Combine plots into a single object and facet_wrap them ----

multi_plot <- cowplot::plot_grid(plots[[1]][[1]], plots[[2]][[1]], plots[[3]][[1]], plots[[4]][[1]],
                                 plots[[1]][[2]], plots[[2]][[2]], plots[[3]][[2]], plots[[4]][[2]],
                                 nrow = 2,
                                 rel_widths = c(1, 1, 1, 1.45), # Relative widths for each column
                                 rel_heights = c(1, 1))

multi_plot

# Run through here (setting plot: 2000*1000 (length*width))




########## Generate separate plots for simulation results ############

# Re-run the code above using the following legend position vector:

# legend_position <- list("right", "right", "right", "right")



















