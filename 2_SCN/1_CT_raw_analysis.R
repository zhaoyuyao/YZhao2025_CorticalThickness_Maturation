## This script is to analyze raw CT
rm(list = ls()) #clears all variables
#  Plots of CT changes of each region is also created
#  Author: Zhao.Yuyao Jul 14 2025
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
setwd("/Users/zhaoyuyao/Desktop/Cohen_lab/EBDS/structural_analysis/1_clean_data") #Sets the working directory to the Desktop

#import csvs
#import csv files
combined_df <- read.csv("CT_train_CS_T1_no_MVM_AAL_globZ2.csv")
matwhole <- combined_df[,4:81]

## plot CT by Yeo2011 (whole brain) ----------------------
index_vis <- c(39:52) 
index_som <- c(1:2,17:18,20,53:54,65:70)
index_lim <- c(5:6,21:22,27:28,37:38,71:72,75:78)
index_ass_pos <- c(7:14,19,29:30,33:34,55:60)
index_ass_neg <- c(3:4,15:16,23:26,31:32,35:36,61:64,73:74)
index_whole <- c(1:78)
mean_vis <- rowMeans(matwhole[, index_vis], na.rm = TRUE)
mean_som <- rowMeans(matwhole[, index_som], na.rm = TRUE)
mean_lim <- rowMeans(matwhole[, index_lim], na.rm = TRUE)
mean_ass_pos <- rowMeans(matwhole[, index_ass_pos], na.rm = TRUE)
mean_ass_neg <- rowMeans(matwhole[, index_ass_neg], na.rm = TRUE)
mean_whole <- rowMeans(matwhole[, index_whole], na.rm = TRUE)
mat_mean <- cbind(combined_df[, 1:3], mean_vis, mean_som, mean_lim, mean_ass_pos,mean_ass_neg,mean_whole)
length(unique(mat_mean$Case))
# Reshape the data from wide to long format for plotting
long_whole <- mat_mean %>%
  pivot_longer(cols = c("mean_vis", "mean_som",  "mean_lim", "mean_ass_pos", "mean_ass_neg", "mean_whole"),
               names_to = "Measure",
               values_to = "Value")
long_whole$Measure <- as.factor(long_whole$Measure)
long_whole$Measure <- factor(long_whole$Measure, levels = c("mean_vis", "mean_som",  "mean_lim", "mean_ass_pos", "mean_ass_neg", "mean_whole"))
levels(long_whole$Measure) # Check the new order

## anova ----------
# Factors & references
long_whole$Sex     <- relevel(factor(long_whole$Sex), ref = "2")
long_whole$Measure <- relevel(factor(long_whole$Measure), ref = "mean_whole")
long_whole$Year    <- factor(long_whole$Year, levels = c("0","1","2","4","6","8","10"))

# (Optional) Type-III setup
options(contrasts = c("contr.sum", "contr.poly"))

# Mixed model (start with intercept-only; consider adding Year slope if warranted)
lmer_whole <- lmer(Value ~ Sex + Measure * Year + (1 | Case), data = long_whole)
summary(lmer_whole)
anova(lmer_whole, type = 3)  # optional

# Pairwise (marginal)
emm_measure <- emmeans(lmer_whole, ~ Measure)
# Pairwise vs global (mean_whole) only
contrast_vs_global <- contrast(emm_measure, method = list(
  "limbic - global"         = c(0, 0,  1, 0, 0, -1),
  "ass_pos - global"        = c(0, 0,  0, 1, 0, -1),
  "ass_neg - global"        = c(0, 0,  0, 0, 1, -1),
  "visual - global"         = c(1, 0,  0, 0, 0, -1),
  "somatomotor - global"    = c(0, 1,  0, 0, 0, -1)
), adjust = "fdr")
means_by_network <- summary(emm_measure) %>%
  select(Measure, emmean, SE, df, lower.CL, upper.CL) %>%
  arrange(Measure)

print(means_by_network)
print(summary(contrast_vs_global))

emm_year <- emmeans(lmer_whole, ~ Year)
pairs(emm_year, adjust = "fdr")

# Simple effects
pairs(emmeans(lmer_whole, ~ Year | Measure),   adjust = "fdr")
pairs(emmeans(lmer_whole, ~ Measure | Year),   adjust = "fdr")



# Additional analysis of yearly changes and amplitude following "1_CT_rawtraj_AAL_Yeo"
## ============================================================
## (B) “Amplitude of yearly changes”: 0–1 vs. later intervals
##     - absolute annualized change in mean_whole
## ============================================================
MEASURE_FOR_CHANGE <- "mean_whole"

# Keep the chosen measure and make Year numeric for interval math
lw_change <- long_whole %>%
  filter(Measure == MEASURE_FOR_CHANGE) %>%
  mutate(Year_num = as.numeric(as.character(Year))) %>%
  arrange(Case, Year_num)

# Compute adjacent-interval annualized absolute changes per subject
changes <- lw_change %>%
  group_by(Case) %>%
  arrange(Year_num, .by_group = TRUE) %>%
  mutate(
    Value_prev = dplyr::lag(Value),
    Year_prev  = dplyr::lag(Year_num),
    gap_years  = Year_num - Year_prev,
    annualized_abs_change = ifelse(!is.na(Value_prev) & gap_years > 0,
                                   abs(Value - Value_prev) / gap_years,
                                   NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(annualized_abs_change))

# Tag early (0–2) vs later (>=2)
changes <- changes %>%
  mutate(period = ifelse(Year_prev == 0 & Year_num == 2, "early_0to2", "later_2to10"))

# Per subject, summarize early and later (mean of available later intervals)
by_subj <- changes %>%
  group_by(Case, period) %>%
  summarize(mean_change = mean(annualized_abs_change, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_change)

# Paired t-test across subjects (subjects with both early_0to2 and later_2to10)
paired_dat <- by_subj %>% filter(!is.na(early_0to2) & !is.na(later_2to10))
t_change <- t.test(paired_dat$early_0to2, paired_dat$later_2to10, paired = TRUE)
t_change
# => Use t_change$statistic and t_change$p.value in your sentence.
# You can also report means:
early_mean  <- mean(paired_dat$early_0to2)
later_mean  <- mean(paired_dat$later_2to10)

## ============================================================
## (C) “Age 1 has the highest individual variance”: test across years
##     We do this on mean_whole; swap Measure if needed.
## ============================================================
lw_global <- long_whole %>%
  filter(Measure == MEASURE_FOR_CHANGE) %>%
  droplevels()

# Variance by year (descriptive)
var_by_year <- lw_global %>%
  group_by(Year) %>%
  summarize(n = dplyr::n(), variance = var(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(variance))
var_by_year
# Identify the max year
year_with_max_var <- var_by_year$Year[1]

# Formal heterogeneity of variance test (robust)
lev <- leveneTest(Value ~ Year, data = lw_global, center = mean)  # F and p-value
lev
# => Report lev[1,"F value"] and lev[1,"Pr(>F)"]

# (Optional) If you want a classic parametric option:
# bartlett.test(Value ~ Year, data = lw_global)

## ============================================================
## (D) Nice printing of the key numbers you'll paste into prose
## ============================================================
cat("\n--- Key numbers for write-up ---\n")
cat(sprintf("Annualized abs change (0-2):  mean = %.4f mm/year\n", early_mean))
cat(sprintf("Annualized abs change (2-10): mean = %.4f mm/year\n", later_mean))
cat(sprintf("Paired t-test early vs later: t = %.3f, df = %d, p = %.4g\n",
            unname(t_change$statistic), unname(t_change$parameter), t_change$p.value))

cat(sprintf("\nVariance by Year (mean_whole):\n"))
print(var_by_year)

cat(sprintf("\nLevene test (Value ~ Year): F = %.3f, df %d,%d, p = %.4g\n",
            lev[1,"F value"], lev[1,"Df"], lev[2,"Df"], lev[1,"Pr(>F)"]))

cat("\nNetwork-level marginal means (mm):\n")
print(means_by_network)

cat("\nContrasts vs global (FDR):\n")
print(summary(contrast_vs_global))
