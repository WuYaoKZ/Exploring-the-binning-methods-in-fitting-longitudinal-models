############ Generate separate and merged plots for coefficients for NHANES application ############

rm(list=ls())
library(tidyverse)
library(dplyr)
library(grid)
library(ggpubr)
library(gridExtra)



# Read in the NHANES application results ----

res_list <- readRDS("./nhanes_res_data/res_list.rds")

# str(res_list[[1]])  # Model fitting results
# str(res_list[[2]])  # Model Fitting Time
# str(res_list[[3]])  # Plots for coef estimates

plot_data <- res_list[[3]]



# Generate separate plots for coefficient estimates for each binning method ----

bin_size <- c(1,10,10,10) # unbin method, and 3 binning methods
method_idx <- 1:4 # 4 methods in total
coef_idx <- 1:10  # 10 coefficients in the model


titles <- c("(Intercept)", "Day of Wear", "Gender (Female)", "Age", 
            "Day of Week (Monday)", "Day of Week (Tuesday)", 
            "Day of Week (Wednesday)", "Day of Week (Thursday)", 
            "Day of Week (Friday)", "Day of Week (Saturday)")  # Note: Sunday is the reference

plots <- list()  # a list to store separate plots

for(i in method_idx){  # loop over 4 methods
  
  plots[[i]] <- list()  # 4 separate lists for 4 methods
  
  bin = bin_size[i]  # bin size
  
  for(j in coef_idx){  # loop over 10 coefficients
    
    plots[[i]][[j]] <- 
      plot_data[[i]][[j]] %>% 
      ggplot() +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_ribbon(aes(x = s, ymax = CI.upper.joint, ymin = CI.lower.joint), fill = "gray20", alpha = 0.2) +
      geom_ribbon(aes(x = s, ymax = CI.upper.pointwise, ymin = CI.lower.pointwise), fill = "gray10", alpha = 0.4) +
      geom_line(aes(x = s, y = beta.hat, color = "Estimate"), alpha = 1, linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_colour_manual(name="", values=c("Estimate"="black")) +
      scale_x_continuous(
        # Setting breaks at every 6 hour
        breaks = seq(0, 1440/bin, by = 360/bin), 
        labels = function(x) {
          hour <- c(0,6,12,18,24)
          minute <- 0
          sprintf("%02d:%02d", hour, minute)  # Format as HH:MM
        }
      )+
      labs(x = "Time of day", y = bquote(paste(beta[.(j)], "(s)")),
           title = titles[j]) +
      theme(legend.position = "none")
    
  }
}

# Check plots for each binning method

# plots[[1]]  # unbin method
# plots[[2]]  # bin and average method
# plots[[3]]  # bin center method
# plots[[4]]  # bin and stack method



# Merge separate coefficient estimates plots between each binning methods ----

title1 = text_grob("Unbin Method", size = 24, face = "bold")
plot1_merged <- grid.arrange(grobs=plots[[1]], nrow=2, top=title1)

title2 = text_grob("Bin and Average Method", size = 24, face = "bold")
plot2_merged <- grid.arrange(grobs=plots[[2]], nrow = 2, top=title2)

title3 = text_grob("Bin Center Method", size = 24, face = "bold")
plot3_merged <- grid.arrange(grobs=plots[[3]], nrow = 2, top=title3)

title4 = text_grob("Bin and Stack Method", size = 24, face = "bold")
plot4_merged <- grid.arrange(grobs=plots[[4]], nrow = 2, top=title4)

### Merge all plots ----

plot_all_merged <- grid.arrange(plot1_merged, plot2_merged, plot3_merged, plot4_merged, ncol = 1)



# Create selected plots (Unbin + Bin and Average) for presentation ----

### Merge Plots
title1 = text_grob("Coefficients Estimates from Unbin Method", size = 16, face = "bold")
plot1 <- grid.arrange(plots[[1]][[1]], plots[[1]][[2]], plots[[1]][[3]], plots[[1]][[4]], plots[[1]][[5]], plots[[1]][[10]],
                      ncol = 2, top=title1)
title2 = text_grob("Coefficients Estimates from Bin and Average Method (Bin Size = 10)", size = 16, face = "bold")
plot2 <- grid.arrange(plots[[2]][[1]], plots[[2]][[2]], plots[[2]][[3]], plots[[2]][[4]], plots[[2]][[5]], plots[[2]][[10]],
                      ncol = 2, top=title2)

# plot1
# plot2









