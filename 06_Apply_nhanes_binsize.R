########## Implement Binning methods using FUI in the NHANES study for different bin sizes ##########

rm(list=ls())
library(tidyverse)
library(refund)
library(dplyr)
library(fastFMM)
library(lme4)
set.seed(1234)

# Read in NHANES multi-level data ----

new_nhanes <- readRDS("./data/nhanes_fda_with_r_ml.rds")
#length(unique(new_nhanes$SEQN))  # Original dataset contains 12610 individuals



# Focus on the young adults ----
# (age between 18 and 30 years old)

select_idx <- new_nhanes$age>=18 & new_nhanes$age<=30
cut_nhanes <- new_nhanes[select_idx,]
cut_nhanes$dayofweek <- factor(cut_nhanes$dayofweek, 
                               levels = c("1","2","3","4","5","6","7"), 
                               labels = c("MON","TUE","WED","THU","FRI","SAT","SUN"))



### Unbin method ----

ptm <- proc.time()
fit1 <- fui(MIMS ~ dayofwear + gender + age + dayofweek + (1|SEQN), 
            data = cut_nhanes, family = "gaussian", analytic = FALSE, num_boots = 50)
time1 <- (proc.time() - ptm)[3]

plot1 <- plot_fui(fit1, num_row = 2, return = TRUE)

# For the next steps, apply four Binning methods to the data and generate a 4 by 4 plot which shows then parameter estimates



# Define a function that loops over 3 binning methods ----

bin_res <- function(bin_size){
  
  L <- ncol(cut_nhanes$MIMS)
  cat("Bin Size: ", bin_size, "\n")
  L_bin <- L/bin_size
  cat("L Bin: ", L_bin, "\n")
  
  ### Bin average method ----
  
  bin_ave_nhanes <- cut_nhanes
  bin_ave_nhanes$bin_ave <- matrix(NA, nrow(bin_ave_nhanes), L_bin)
  
  ptm <- proc.time()
  for(j in 1:L_bin){
    bin_ave_nhanes$bin_ave[,j] <- rowMeans(bin_ave_nhanes$MIMS[,(bin_size*(j-1)+1):(bin_size*j)])
  }
  time2.1 <- (proc.time() - ptm)[3]
  
  print("Model: bin and average")    
  ptm <- proc.time()
  fit2 <- fui(bin_ave ~ dayofwear + gender + age + dayofweek + (1|SEQN), 
              data = bin_ave_nhanes, family = "gaussian")
  time2.2 <- (proc.time() - ptm)[3]
  time2 <- c(time2.1 + time2.2, time2.1, time2.2)
  
  plot2 <- plot_fui(fit2, num_row = 2, return = TRUE)
  
  
  ### Bin center method ----
  
  bin_center_nhanes <- cut_nhanes
  bin_center_nhanes$bin_center <- matrix(NA, nrow(bin_center_nhanes), L_bin)
  
  ptm <- proc.time()
  for(i in 1:nrow(bin_center_nhanes)){
    bin_center_nhanes$bin_center[i,] <- bin_center_nhanes$MIMS[i, seq(1:L_bin)*bin_size - floor(bin_size/2)]
  }
  time3.1 <- (proc.time() - ptm)[3]
  
  print("Model: bin center")  
  ptm <- proc.time()
  fit3 <- fui(bin_center ~ dayofwear + gender + age + dayofweek + (1|SEQN), 
              data = bin_center_nhanes, family = "gaussian")
  time3.2 <- (proc.time() - ptm)[3]
  time3 <- c(time3.1 + time3.2, time3.1, time3.2)
  
  plot3 <- plot_fui(fit3, num_row = 2, return = TRUE)
  
  
  
  ### Bin and stack method ----
  
  New_SEQN <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
  New_dayofwear <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
  New_dayofweek <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
  New_gender <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
  New_age <- matrix(NA, nrow(cut_nhanes)*bin_size, 1)
  
  New_SEQN <- as.factor(rep(cut_nhanes$SEQN, each = bin_size))
  New_dayofwear <- rep(cut_nhanes$dayofwear, each = bin_size)
  New_dayofweek <- rep(cut_nhanes$dayofweek, each = bin_size)
  New_gender <- rep(cut_nhanes$gender, each = bin_size)
  New_age <- rep(cut_nhanes$age, each = bin_size)
  
  New_MIMS <- matrix(NA, nrow(cut_nhanes)*bin_size, L_bin)
  
  ptm <- proc.time()
  for(i in 1:nrow(cut_nhanes)){
    for(j in 1:L_bin){
      New_MIMS[(bin_size*(i-1)+1):(bin_size*i), j] <- cut_nhanes$MIMS[i, (bin_size*(j-1)+1):(bin_size*j)]
    }
  }
  bin_stack_nhanes <- data.frame(SEQN=New_SEQN, dayofwear=New_dayofwear, dayofweek=New_dayofweek, 
                                 gender=New_gender, age=New_age, bin_stack=I(New_MIMS))
  time4.1 <- (proc.time() - ptm)[3]
  
  print("Model: bin and stack")
  ptm <- proc.time()
  fit4 <- fui(bin_stack ~ dayofwear + gender + age + dayofweek + (1|SEQN), 
              data = bin_stack_nhanes, family = "gaussian")
  time4.2 <- (proc.time() - ptm)[3]
  time4 <- c(time4.1 + time4.2, time4.1, time4.2)
  
  plot4 <- plot_fui(fit4, num_row = 2, return = TRUE)
  
  # Save and return the results
  return(
    list(bin_size = bin_size, 
         fit = list(fit2,fit3,fit4),      # model fitting results
         time = list(time2,time3,time4),  # computing time
         plot = list(plot2,plot3,plot4))  # coefficient estimates plots
  )
  
}



# Get model fitting results for different bin sizes ----

bin_60 <- bin_res(bin_size = 60)
bin_50 <- bin_res(bin_size = 50)  # the number of bins is not integer
bin_40 <- bin_res(bin_size = 40)
bin_30 <- bin_res(bin_size = 30)
bin_20 <- bin_res(bin_size = 20)
bin_10 <- bin_res(bin_size = 10)

bin_res_list <- list(bin_10, bin_20, bin_30, bin_40, bin_50, bin_60)  # results for all binning methods

unbin <- list(bin_size = 1, fit = fit1, time = time1, plot = plot1)  # results for unbin method

all_res_list <- list(unbin = unbin, 
                     bin_10 = bin_10, bin_20 = bin_20, 
                     bin_30 = bin_30, bin_40 = bin_40, 
                     bin_50 = bin_50, bin_60 = bin_60)  # results for unbin and binning methods

# Save results

saveRDS(bin_res_list, "./nhanes_res_data/bin_res_list.rds")
saveRDS(all_res_list, "./nhanes_res_data/all_res_list.rds")



# Focus on "bin and average" method and small bin sizes ----

bin_ave_res <- function(bin_size){
  
  # Reminder of bin size and number of bins
  L <- ncol(cut_nhanes$MIMS)
  cat("Bin Size: ", bin_size, "\n")
  L_bin <- L/bin_size
  cat("L Bin: ", L_bin, "\n")
  
  # Bin average method
  bin_ave_nhanes <- cut_nhanes
  bin_ave_nhanes$bin_ave <- matrix(NA, nrow(bin_ave_nhanes), L_bin)
  
  ptm <- proc.time()
  for(j in 1:L_bin){
    bin_ave_nhanes$bin_ave[,j] <- rowMeans(bin_ave_nhanes$MIMS[,(bin_size*(j-1)+1):(bin_size*j)])
  }
  
  print("Model: bin and average")    
  fit2 <- fui(bin_ave ~ dayofwear + gender + age + dayofweek + (1|SEQN), 
              data = bin_ave_nhanes, family = "gaussian")
  time2 <- (proc.time() - ptm)[3]
  
  plot2 <- plot_fui(fit2, num_row = 2, return = TRUE)
  
  # Save the results
  return(
    list(bin_size = bin_size, 
         fit = fit2,
         time = time2,
         plot = plot2)
  )
  
}

# Exploring selected bin sizes (3,5,8,10,20,30,60)

bin_ave_3 <- bin_ave_res(bin_size = 3)
bin_ave_5 <- bin_ave_res(bin_size = 5)
bin_ave_8 <- bin_ave_res(bin_size = 8)
bin_ave_10 <- bin_ave_res(bin_size = 10)
bin_ave_20 <- bin_ave_res(bin_size = 20)
bin_ave_30 <- bin_ave_res(bin_size = 30)
bin_ave_60 <- bin_ave_res(bin_size = 60)

bin_ave_3$time   # 2.51 hours = 150.7 minutes
bin_ave_5$time   # 0.92 hours = 55.2 minutes
bin_ave_8$time   # 0.37 hours = 22.4 minutes
bin_ave_10$time  # 0.24 hours = 14.4 minutes
bin_ave_20$time  # 0.07 hours = 4.0 minutes
bin_ave_30$time  # 0.03 hours = 1.9 minutes
bin_ave_60$time  # 0.01 hours = 0.5 minute

# Save the results
saveRDS(bin_ave_3, "./nhanes_res_data/bin_ave_3.rds")
saveRDS(bin_ave_5, "./nhanes_res_data/bin_ave_5.rds")
saveRDS(bin_ave_8, "./nhanes_res_data/bin_ave_8.rds")
saveRDS(bin_ave_10, "./nhanes_res_data/bin_ave_10.rds")
saveRDS(bin_ave_20, "./nhanes_res_data/bin_ave_20.rds")
saveRDS(bin_ave_30, "./nhanes_res_data/bin_ave_30.rds")
saveRDS(bin_ave_60, "./nhanes_res_data/bin_ave_60.rds")



# Plots of estimated coefficients and computing time for "bin and average" method ----

rm(list=ls())
library(grid)
library(ggpubr)
library(gridExtra)



### Read in model fitting results data ----

all_res_list <- readRDS("./nhanes_res_data/all_res_list.RDS")
unbin_res_list <- all_res_list[[1]]

bin_ave_3 <- readRDS("./nhanes_res_data/bin_ave_3.RDS")
bin_ave_5 <- readRDS("./nhanes_res_data/bin_ave_5.RDS")
bin_ave_8 <- readRDS("./nhanes_res_data/bin_ave_8.RDS")
bin_ave_10 <- readRDS("./nhanes_res_data/bin_ave_10.RDS")
bin_ave_20 <- readRDS("./nhanes_res_data/bin_ave_20.RDS")
bin_ave_30 <- readRDS("./nhanes_res_data/bin_ave_30.RDS")
bin_ave_60 <- readRDS("./nhanes_res_data/bin_ave_60.RDS")

plot_unbin <- unbin_res_list$plot
plot_bin_3 <- bin_ave_3$plot
plot_bin_5 <- bin_ave_5$plot
plot_bin_8 <- bin_ave_8$plot
plot_bin_10 <- bin_ave_10$plot
plot_bin_20 <- bin_ave_20$plot
plot_bin_30 <- bin_ave_30$plot
plot_bin_60 <- bin_ave_60$plot

plot_data <- list(plot_unbin, plot_bin_3, plot_bin_5, plot_bin_8, 
                  plot_bin_10, plot_bin_20, plot_bin_30, plot_bin_60)



### Plot all coefficient estimates for different bin sizes ----

bin_size <- c(1,3,5,8,10,20,30,60)
binsize_idx <- 1:length(bin_size)  # 8 binsizes
titles <- c("(Intercept)", "Day of Wear", "Gender (Female)", "Age",
            "Day of Week (Monday)", "Day of Week (Tuesday)",
            "Day of Week (Wednesday)", "Day of Week (Thursday)",
            "Day of Week (Friday)", "Day of Week (Saturday)")
coef_idx <- 1:length(titles)     # 10 Parameter Coefficients

plots <- list()

for(i in binsize_idx){  # loop over 4 methods
  
  plots[[i]] <- list()
  
  bin = bin_size[i]
  
  for(j in coef_idx){  # loop over 10 coeffs
    
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
        breaks = seq(0, 1440/bin, by = 360/bin),  # Setting breaks at every 6 hour
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

title1 = text_grob("Unbin Method", size = 24, face = "bold")
plot1_merged <- grid.arrange(grobs=plots[[1]], nrow=2, top=title1)

title2 = text_grob("Bin and Average (Bin=3)", size = 24, face = "bold")
plot2_merged <- grid.arrange(grobs=plots[[2]], nrow = 2, top=title2)

title3 = text_grob("Bin and Average (Bin=5)", size = 24, face = "bold")
plot3_merged <- grid.arrange(grobs=plots[[3]], nrow = 2, top=title3)

title4 = text_grob("Bin and Average (Bin=8)", size = 24, face = "bold")
plot4_merged <- grid.arrange(grobs=plots[[4]], nrow = 2, top=title4)

title5 = text_grob("Bin and Average (Bin=10)", size = 24, face = "bold")
plot5_merged <- grid.arrange(grobs=plots[[5]], nrow = 2, top=title5)

title6 = text_grob("Bin and Average (Bin=20)", size = 24, face = "bold")
plot6_merged <- grid.arrange(grobs=plots[[6]], nrow = 2, top=title6)

title7 = text_grob("Bin and Average (Bin=30)", size = 24, face = "bold")
plot7_merged <- grid.arrange(grobs=plots[[7]], nrow=2, top=title7)

title8 = text_grob("Bin and Average (Bin=60)", size = 24, face = "bold")
plot8_merged <- grid.arrange(grobs=plots[[8]], nrow = 2, top=title8)

# Merge all plots

plot_all_merged <- grid.arrange(plot1_merged, plot2_merged, plot3_merged, plot4_merged,
                                plot5_merged, plot6_merged, plot7_merged, plot8_merged, 
                                ncol = 1)



### Plot coefficient estimates for a single covariate (e.g. "Day of wear") ----

# new titles indicating different bin sizes
titles <- c("Unbin", "Bin Size = 3", 
            "Bin Size = 5", "Bin Size = 8",
            "Bin Size = 10", "Bin Size = 20",
            "Bin Size = 30", "Bin Size = 60")

plots <- list()

for(i in binsize_idx){  # loop over 8 bin sizes
  
  plots[[i]] <- list()
  
  bin = bin_size[i]
  
  for(j in coef_idx){  # loop over 10 coeffs
    
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
        breaks = seq(0, 1440/bin, by = 360/bin),  # Setting breaks at every 6 hour
        labels = function(x) {
          hour <- c(0,6,12,18,24)
          minute <- 0
          sprintf("%02d:%02d", hour, minute)  # Format as HH:MM
        }
      )+
      labs(x = "Time of day", y = bquote(paste(beta[.(j)], "(s)")),
           title = titles[i]) +
      theme(legend.position = "none")
    
  }
}

# use the 2nd coefficient "day of wear" for illustration

# plots[[1]][[2]]  # Unbin
# plots[[2]][[2]]  # Bin=3
# plots[[3]][[2]]  # Bin=5
# plots[[4]][[2]]  # Bin=8
# plots[[5]][[2]]  # Bin=10
# plots[[6]][[2]]  # Bin=20
# plots[[7]][[2]]  # Bin=30
# plots[[8]][[2]]  # Bin=60

title1 = text_grob('Coefficients Estimates for "Day of Wear"', size = 16, face = "bold")
dayofwear_bybinsize <- grid.arrange(plots[[1]][[2]], plots[[2]][[2]],
                                    plots[[3]][[2]], plots[[4]][[2]],
                                    plots[[5]][[2]], plots[[6]][[2]],
                                    plots[[7]][[2]], plots[[8]][[2]],
                                    ncol = 3, top=title1)



### Compare model fitting time for different bin sizes ----

t1 <- unbin_res_list$time     # 6.3 hours
t2 <- bin_ave_3$time          # 2.50 hours
t3 <- bin_ave_5$time          # 0.91 hours
t4 <- bin_ave_8$time          # 0.36 hours
t5 <- bin_ave_10$time[[1]]    # 0.24 hours
t6 <- bin_ave_20$time[[1]]    # 0.07 hours
t7 <- bin_ave_30$time[[1]]    # 0.03 hours
t8 <- bin_ave_60$time[[1]]    # 0.01 hours

time_table <- data.frame(time = c(t1,t2,t3,t4,t5,t6,t7,t8), row.names = c("1", "3", "5", "8",
                                                                          "10", "20","30", "60"))
ref_time <- time_table[1, 1]  # reference computing time: unbin method
time_table$percentage <- (time_table$time / ref_time) * 100  # scale the number

# Plots of computing time for different bin sizes

time_table_plot <- cbind(binsize = c(1,3,5,8,10,20,30,60), time_table)

time_table_plot %>% 
  ggplot(aes(x=binsize, y=percentage)) +
  geom_point(col="blue") +
  geom_line() +
  geom_hline(yintercept=10, col="red", linetype=2) +
  xlab("Bin Size (minute)") +
  ylab("Percentage (%)") +
  labs(title = "Percentages of Computing Time for different Bin Sizes") +
  scale_x_continuous(breaks = c(1,3,5,8,10,20,30,60)) +
  scale_y_continuous(breaks = c(0,10,25,50,75,100)) +
  theme_bw()

# Table for the computing time for different bin sizes

colnames(time_table) <- c("Time (second)", "Percentages (%)")
time_table <- round(time_table, digits = 1)

library(magrittr)
library(flextable)
time_table %>% 
  tibble::rownames_to_column() %>% 
  flextable() %>%
  set_header_labels(rowname = "Bin Sizes") %>% 
  theme_zebra() %>% autofit()









