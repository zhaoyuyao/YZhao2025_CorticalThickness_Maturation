## This script is to reorganize data (as used for SCN) for following sbMCN analysis
#  Author: Zhao.Yuyao Mar 10 2025

rm(list = ls()) #clears all variables
# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("...") 

#import csv files
combined_df <- read.csv("CT_test_X_XCT_T2_no_MVM_AAL_thresholdZ2.csv")

#Whole Brain
index <- c(1,3,2,4:81) 
matReorder1 <- combined_df[index]
table(matReorder1$Year)
unique_sub<-matReorder1 %>% group_by(Case, Sex) %>% slice_sample(n = 1)
table(unique_sub$Sex)
length(unique_sub$Case)
write.csv(matReorder1, file = ".../sbMCN_avg_z_sort_WB_no_MVM_Y0to10.csv")

