## Modified for EBDS CT SCN analysis
## Outlier removal based on overall mean CT
## Yuyao Zhao, Aug 21 2024; updated Aug 14 2025

rm(list = ls())

library(dplyr)
library(ggplot2)

setwd("...")
save.dir <- ".../3_WB_matprep_regoverallave_globZ2_no_MVM_AAL"
dir.create(save.dir)
dir.create(fig.dir)

file_name <- "CT_only_CS_T1_no_MVM_AAL"
data <- read.csv("CT_only_CS_T1_no_MVM_AAL.csv")

# Convert Sex to numeric codes
data$Sex[data$Sex == "M"] <- "1"
data$Sex[data$Sex == "F"] <- "0"
data$Sex <- as.numeric(data$Sex)

# List of years to process
years <- c(0, 1, 2, 4, 6, 8, 10)

# Store outputs
outliers_all <- list()
resid_list <- list()
matID_list <- list()

for (yr in years) {
  
  # Filter year
  my_data <- data %>% filter(Year == yr)
  
  # CT matrix
  allCT <- my_data[, 4:81]
  
  # Overall mean CT z-score
  my_data$CTmeans <- rowMeans(allCT, na.rm = TRUE)
  my_data$CTmeans_z <- scale(my_data$CTmeans)
  CTmeanscol <- colMeans(allCT, na.rm = TRUE) #calculates mean for each ROI column
  my_data$OverallMT <-  mean(CTmeanscol) #calculates overall mean CT and recodes column with mean
  
  
  overall_outlier_mask <- abs(my_data$CTmeans_z) >= 2
  combined_mask <- overall_outlier_mask
  
  # Save outliers for checking
  if (any(combined_mask)) {
    outliers_all[[paste0("T", yr)]] <- my_data[combined_mask, ]
  }
  
  # Keep only non-outliers
  my_data <- my_data[!combined_mask, ]
  allCT <- my_data[, 4:81]
  
  # Regression: remove effects of sex + mean overall CT from each ROI
  reg <- lapply(4:81, function(i) {
    resid(lm(my_data[, i] ~ my_data$Sex + my_data$OverallMT))
  })
  
  # Combine residuals into a matrix
  resid_all <- do.call(cbind, reg)
  resid_list[[paste0("T", yr)]] <- resid_all
  matID_list[[paste0("T", yr)]] <- data.frame(Case = my_data$Case, resid_all)
  
  # Save outputs
  write.csv(resid_all, file = file.path(save.dir, paste0("T", which(years == yr), "_DK_resid_age.csv")), row.names = FALSE)
  write.csv(matID_list[[paste0("T", yr)]], file = file.path(save.dir, paste0("T", which(years == yr), "_DK_ID_resid_age.csv")), row.names = FALSE)
  
  # Save correlation matrix
  write.csv(cor(resid_all), file = file.path(save.dir, paste0("T", which(years == yr), "_WB_reg_gender_WBCT.csv")))

# Combine all non-outlier data
combined_df <- bind_rows(lapply(names(resid_list), function(nm) {
  idx <- which(names(resid_list) == nm)
  yr <- years[idx]
  my_data <- data %>% filter(Year == yr)
  
  # Remove outliers again for combined_df
  outlier_cases <- outliers_all[[paste0("T", yr)]]$Case
  my_data <- my_data[!(my_data$Case %in% outlier_cases), ]
  my_data
}))

}
write.csv(combined_df, file = paste0(file_name, "regoverallave_globZ2_train.csv"), row.names = FALSE)
# later renamed to "CT_train_CS_T1_no_MVM_AAL_globZ2"
# Save outliers to check later 
outliers_df <- bind_rows(outliers_all, .id = "Timepoint")
write.csv(outliers_df, file = paste0(file_name, "regoverallave_globZ2_train_outliers.csv"), row.names = FALSE)
# later renamed to "CT_train_CS_T1_no_MVM_AAL_globZ2_outliers"