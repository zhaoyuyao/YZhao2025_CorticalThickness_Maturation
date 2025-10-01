## this code is to merge and standardize all the behavioral variebles into one file


# Clear the environment
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(corrr)
# Load files
behav_data1 <- read.csv(".../Flanker_10older_no_MVM.csv")
behav_data2 <- read.csv(".../DCCS_10older_no_MVM.csv")
behav_data3 <- read.csv(".../LS_10older_no_MVM.csv")
behav_data4 <- read.csv(".../PS_10older_no_MVM.csv")
behav_data5 <- read.csv(".../SSP_10older_no_MVM.csv")

## ========= 1) Merge all behavioral data and select oldest timepoint =========

# Merge the behavioral data frames
behav_data <- merge(merge(merge(merge(behav_data1, behav_data2,
                                      by = c("Case","SubjectID","Sex","Year"),all=TRUE),
                                behav_data3, by = c("Case","SubjectID","Sex","Year"),all=TRUE),
                          behav_data4, by = c("Case","SubjectID","Sex","Year"),all=TRUE),
                    behav_data5, by = c("Case","SubjectID","Sex","Year"),all=TRUE)

# Select the oldest timepoint for each SubjectID
behav_data_oldest <- behav_data[which(behav_data$Year>10&behav_data$Year<18),] %>%
  group_by(SubjectID) %>%
  filter(Year == max(Year, na.rm=TRUE)) %>%
  ungroup()



## ========= 2) Standardize variebles ===========

# Inspect
length(unique(behav_data_oldest$SubjectID))
table(behav_data_oldest$Year)
table(behav_data_oldest$Sex)
head(behav_data_oldest)
str(behav_data_oldest)

# Scale your behavioral variables
behav_data_oldest$ssp_length_scaled <- scale(behav_data_oldest$ssp_length)
behav_data_oldest$nt_fic_as_scaled  <- scale(behav_data_oldest$nt_fic_as)
behav_data_oldest$nt_dccs_as_scaled <- scale(behav_data_oldest$nt_dccs_as)
behav_data_oldest$nt_ps_as_scaled   <- scale(behav_data_oldest$nt_ps_as)
behav_data_oldest$nt_ls_as_scaled   <- scale(behav_data_oldest$nt_ls_as)

# Composite Score using PCA
complete_idx <- complete.cases(behav_data_oldest[, c("nt_ls_as_scaled","ssp_length_scaled", "nt_fic_as_scaled", "nt_dccs_as_scaled")])
# Initialize composite score as NA
behav_data_oldest$composite_score <- NA
pca_result <- prcomp(behav_data_oldest[complete_idx, c("nt_ls_as_scaled","ssp_length_scaled", "nt_fic_as_scaled", "nt_dccs_as_scaled")], scale. = TRUE)
behav_data_oldest$composite_score[complete_idx] <- -pca_result$x[, 1]  # Use the first principal component
cor(behav_data_oldest[complete_idx, c("composite_score","ssp_length_scaled", "nt_fic_as_scaled", "nt_dccs_as_scaled","nt_ps_as_scaled","nt_ls_as_scaled")], use="complete.obs")
summary(pca_result)
pca_result$rotation[,1]

write.csv(behav_data_oldest, file = ".../behav_10oldest_train_merged_scaled_no_MVM.csv")
