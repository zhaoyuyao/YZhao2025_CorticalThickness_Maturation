#------------------------------
# Project: EBDS
# Title: SCN metrics overall ANOVA and t-test between networks
# Author: Yuyao Zhao
# Date: Sep 23 2025
#------------------------------
rm(list = ls())

# --- Paths ---
setwd('...')
# Load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

library(tidyverse)
library(emmeans)

# Load aggregated data
df <- read_csv("SCN_metrics_long_globZ2_train.csv") %>%
  filter(Measure %in% c("Seg", "GE", "Degree")) %>%
  filter(Network != "WB") %>%   # <-- drop whole-brain
  mutate(
    Year = factor(Year),
    Network = factor(Network),
    # make sure Measure is a LABEL, not an index
    Measure = as.character(Measure)
  )

# ---- Function: ANOVA + emmeans ----
run_anova_emm <- function(data, measure_name) {
  cat("\n===== ", measure_name, " =====\n")
  
  # Factorial ANOVA (Year × Network)
  mod <- aov(Mean ~ Network, data = data)
  print(summary(mod))
  
  # Estimated marginal means for Network within Year
  emm <- emmeans(mod, ~ Network)
  print(emm)
  
  # Pairwise comparisons among Networks (within each Year)
  cmp <- pairs(emm, adjust = "fdr")
  print(cmp)
}

# ---- Apply to each measure ----
df %>%
  group_split(Measure) %>%
  walk(~ run_anova_emm(.x, unique(.x$Measure)))
