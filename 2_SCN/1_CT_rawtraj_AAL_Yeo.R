## This script is to plot raw data (as used for SCN) 
rm(list = ls()) #clears all variables
#  Plots of CT changes of each region is also created
#  Author: Zhao.Yuyao Jul 14 2025
library(ggplot2)
library(tidyr)
library(dplyr)
library(emmeans)
setwd("...") #Sets the working directory to the Desktop

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

# Reshape the data from wide to long format for plotting
long_whole <- mat_mean %>%
  pivot_longer(cols = c("mean_vis", "mean_som",  "mean_lim", "mean_ass_pos", "mean_ass_neg", "mean_whole"),
               names_to = "Measure",
               values_to = "Value")
long_whole$Measure <- as.factor(long_whole$Measure)
long_whole$Measure <- factor(long_whole$Measure, levels = c("mean_vis", "mean_som",  "mean_lim", "mean_ass_pos", "mean_ass_neg", "mean_whole"))
levels(long_whole$Measure) # Check the new order
# figure information set up
fig.dpi <- 300; fig.wid <- 13; fig.hei <- 10; fig.fmt <- "png" 
fig.savedir <- getwd()

plot.wholetraj <- 
  long_whole %>% 
  ggplot(aes(x = Year, y = Value, color = Measure, fill = Measure, group = Measure)) +
  geom_point(aes(fill = Measure), position = position_jitter(w = 0, h = 0), size = 0.5, shape = 21) + 
  scale_fill_manual(
    values = c('dodgerblue', 'chartreuse4', 'coral', 'firebrick3', 'dodgerblue3', 'gray'),
    labels = c('Vis', 'Som', 'Lim','Pos', 'Neg', 'WB')
  ) + 
  scale_color_manual(
    values = c('dodgerblue', 'chartreuse4', 'coral', 'firebrick3', 'dodgerblue3', 'gray'),
    labels = c('Vis', 'Som', 'Lim','Pos', 'Neg', 'WB')
  ) +
  scale_y_continuous(breaks = c(2,3,4,5,6), limits =c(1.7,5.3) )+
  scale_x_continuous(labels = as.numeric(long_whole$Year), breaks = long_whole$Year) + 
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.8),
        axis.line.y = element_line(colour = "black", size = 0.8),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15, colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "right"  # Ensure legend is displayed
  ) +
  xlab('Age (years)') + 
  ylab("Cortical Thickness") +
  stat_summary(aes(Year, Value), geom = "point", fun = mean) +  #Simple line plot
  stat_summary(aes(Year, Value), geom = "line",  fun = mean) +
  #  geom_smooth(aes(color = Measure), method = "loess", size = 1.5) + 
  guides(color = guide_legend(title = "Networks"),
         fill = guide_legend(title = "Networks"))
plot.wholetraj

# figure name
fig.name <- paste("plot_whole_traj", ".", fig.fmt, sep = "")
## save figure
ggsave(fig.name, path = fig.savedir, plot.wholetraj, width = fig.wid, height = fig.hei, 
       units = "cm", dpi = fig.dpi, bg = "white")


## anova ----------
long_whole$Measure <- as.factor(long_whole$Measure)
levels(long_whole$Measure)
long_whole$Year <- as.factor(long_whole$Year)
levels(long_whole$Year)
# Fit the model with random intercept for Case (subject)
lmer_whole <- lmer(Value ~ Sex + Measure * Year + (1 | Case), data = long_whole)

# Post-hoc pairwise comparisons
emm_net <- emmeans(lmer_whole, ~ Measure)
pairwise_contrasts <- contrast(emm_net, method = "pairwise", adjust = "fdr")
print(pairwise_contrasts)

emm_year <- emmeans(lmer_whole, ~ Year)
pairwise_contrasts <- contrast(emm_year, method = "pairwise", adjust = "fdr")
print(pairwise_contrasts)

# Post-hoc test simple effects
emm <- emmeans(lmer_whole, ~ Year|Measure )
contrast(emm, method = "pairwise", adjust = "fdr")  # or "holm" or "fdr"

emm <- emmeans(lmer_whole, ~ Measure|Year)
contrast(emm, method = "pairwise", adjust = "fdr")  # or "holm" or "fdr"


