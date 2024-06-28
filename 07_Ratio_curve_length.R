########## Explore the bin size choice based on ratio of length of coefficient estimates ##########

rm(list=ls())
library(tidyverse)
library(dplyr)
library(grid)
library(ggpubr)
set.seed(1234)

# Read in the binning methods results data for different bin sizes ----

all_res_list <- readRDS("./nhanes_res_data/all_res_list.RDS")
unbin_res_list <- all_res_list[[1]]

bin_ave_3 <- readRDS("./nhanes_res_data/bin_ave_3.RDS")
bin_ave_4 <- readRDS("./nhanes_res_data/bin_ave_4.RDS")
bin_ave_5 <- readRDS("./nhanes_res_data/bin_ave_5.RDS")
bin_ave_6 <- readRDS("./nhanes_res_data/bin_ave_6.RDS")
bin_ave_8 <- readRDS("./nhanes_res_data/bin_ave_8.RDS")
bin_ave_10 <- readRDS("./nhanes_res_data/bin_ave_10.RDS")
bin_ave_20 <- readRDS("./nhanes_res_data/bin_ave_20.RDS")
bin_ave_30 <- readRDS("./nhanes_res_data/bin_ave_30.RDS")
bin_ave_60 <- readRDS("./nhanes_res_data/bin_ave_60.RDS")

betaHat_unbin <- unbin_res_list$fit$betaHat
betaHat_bin_3 <- bin_ave_3$fit$betaHat
betaHat_bin_4 <- bin_ave_4$fit$betaHat
betaHat_bin_5 <- bin_ave_5$fit$betaHat
betaHat_bin_6 <- bin_ave_6$fit$betaHat
betaHat_bin_8 <- bin_ave_8$fit$betaHat
betaHat_bin_10 <- bin_ave_10$fit$betaHat
betaHat_bin_20 <- bin_ave_20$fit$betaHat
betaHat_bin_30 <- bin_ave_30$fit$betaHat
betaHat_bin_60 <- bin_ave_60$fit$betaHat



# Calculate the curve length of coefficient estimates for different bin sizes ----

binsize_list <- c(1,3,4,5,6,8,10,20,30,60)
curve_length_bin_list <- array(0, dim = c(length(binsize_list),
                                          dim(betaHat_bin_10)[1]))
betaHat_list <- list(betaHat_unbin,
                     betaHat_bin_3,
                     betaHat_bin_4,
                     betaHat_bin_5,
                     betaHat_bin_6,
                     betaHat_bin_8,
                     betaHat_bin_10,
                     betaHat_bin_20,
                     betaHat_bin_30,
                     betaHat_bin_60)

for(k in 1:length(binsize_list)){  #loop over different bin sizes
  
  bin_size <- binsize_list[k]
  
  for(i in 1:dim(curve_length_bin_list)[2]){  #loop over parameter coefficients
    
    coefs <- as.numeric(betaHat_list[[k]][i,])
    curve_length_bin_list[k,i] <- bin_size/1440 # assume first piece of curve to be horizontal
    
    for (j in 2:length(coefs)){  #loop over points on the functional domain
      
      curve_length_bin_list[k,i] <- curve_length_bin_list[k,i] + sqrt((coefs[j]-coefs[j-1])^2 + (bin_size/1440)^2)  # Normalize the function domain to be 0 to 1
    }
  }
}

# Calculate the ratio of curve lengths between "bin and average" method and "unbin" method ----

ratios <- array(NA, dim = c(length(binsize_list)-1, dim(betaHat_bin_10)[1]))

for(k in 2:dim(curve_length_bin_list)[1]){
  ratios[k-1,] <- curve_length_bin_list[k,] / curve_length_bin_list[1,]
}

rownames(ratios) <- c("Binsize=3", "Binsize=4", "Binsize=5", "Binsize=6", 
                      "Binsize=8", "Binsize=10", "Binsize=20", "Binsize=30", "Binsize=60")
colnames(ratios) <- c("(Intercept)", "Day of Wear", "Female", "Age", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
ratios <- round(ratios, 3)
ratios



# Plot the average curve length ratios for different bin sizes ----

ave_ratios <- rowMeans(ratios) # averaged value for all covariates
ave_ratios <- cbind(ave_ratios, binsize = c(3,4,5,6,8,10,20,30,60))
ave_ratios <- as.data.frame(ave_ratios)

ave_ratio_plot <- 
  ggplot(data=ave_ratios, aes(x=binsize, y=ave_ratios)) +
  geom_line(linewidth=0.8) +
  labs(title = "Average Curve Ratio Over All Coefficient Estimates") +
  xlab("Bin Size (minute)") +
  ylab("Curve Length Ratios") +
  scale_x_continuous(
    breaks = c(5,10,20,30,60)
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_bw()



# Generate separate and merged plots for covariates ----

ratios <- cbind(ratios, binsize = c(3,4,5,6,8,10,20,30,60))
ratios <- as.data.frame(ratios)  # ratio by bin size across all covariates


var_list <- colnames(ratios[,-11])
plot_list <- list()

for(i in 1:length(var_list)){
  
  plot_list[[i]] <- 
    eval(substitute(
      ggplot(data=ratios, aes(x=binsize, y=ratios[,i])) +
        geom_line(linewidth=0.8) +
        labs(title = var_list[i]) +
        xlab("Bin Size (min)") +
        ylab("Curve Length Ratios") +
        scale_x_continuous(
          breaks = c(3,4,5,6,8,10,20,30,60),
          labels = c("3","4","5","6","8","10","20","30","60")  # use actual bin sizes as x-axis ticks
          #,
          #labels = c("1/480","1/240","1/144","1/72","1/48","1/24"),   # change the x-axis ticks to relative measure
          #guide = guide_axis(angle = 45)
          
        ) +
        scale_y_continuous(limits = c(0, 0.5)) +
        theme_bw(), list(i=i)))
  
}

plot_list

# Arrange plots into a 5x2 grid
grid.arrange(grobs = plot_list, ncol = 2)

# Arrange plots into a 3x2 grid (one weekday + one weekend)
sub_plot_list <- list(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                      plot_list[[4]], plot_list[[5]], plot_list[[10]])
title1 <- text_grob("Curve Length Ratio between Bin and Unbin Methods", size = 16, face = "bold")
grid.arrange(grobs = sub_plot_list, ncol = 2, top=title1)


