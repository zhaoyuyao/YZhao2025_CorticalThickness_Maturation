#------------------------------
# Project: EBDS
# Title: Subject distribution figure for EBDS behav data
# Author: Yuyao Zhao
# Date: Jul 22 2025
#------------------------------

rm(list = ls())

# Change directory path
setwd('...')
cur.dir <- getwd()
# Define figure folder
fig.fold <- "R_figrues"
# Combine the current directory path with the new folder name
fig.savedir <- file.path(cur.dir, fig.fold)
# Create the new directory
#dir.create(fig.savedir)

# Loading packages 
packages =  c("ggplot2",'ggpubr','dplyr','tidyverse')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), repos = "http://cran.us.r-project.org")
}
invisible(lapply(packages, library, character.only = TRUE))

# Loading files
df <- read.csv("behav_10oldest_train_merged_scaled_no_MVM.csv")


unique_sub<-df %>% group_by(Case, Sex) %>% slice_sample(n = 1)
table(unique_sub$Sex)
length(unique_sub$Case)
table(df$Year)



# Create a long format dataset
df_long <- df %>%
  pivot_longer(
    cols = c(nt_fic_as, nt_dccs_as, nt_ls_as, nt_ps_as, ssp_length),
    names_to = "Task",
    values_to = "Score"
  ) %>%
  mutate(Task = recode(Task,
                       ssp_length = "SSP",
                       nt_ls_as = "LS",
                       nt_fic_as = "Flanker",
                       nt_dccs_as = "DCCS",
                       nt_ps_as = "PCPS"
                      ))

# Count number of valid scores per Year, Task, and Sex
df_count <- df_long %>%
  filter(!is.na(Score)) %>%
  group_by(Year, Task, Sex) %>%
  summarise(N = n(), .groups = "drop")
# Reverse stacking order: F on bottom, M on top
df_count$Sex <- factor(df_count$Sex, levels = c("M","F"),labels = c("male","female"))
df_count$Task <- factor(df_count$Task,
                        levels = c("SSP","LS", "Flanker", "DCCS", "PCPS"))

# Plot


Sub.info<-ggplot(df_count, aes(x = factor(Year), y = N, fill = Sex, color = Sex)) +
  geom_bar(stat = "identity", width = 0.5) +  # stacked bars by default
  facet_wrap(~ Task, nrow = 1) +  # one panel per task
  labs(x = "Age (Years)", y = "#Subject") +
  scale_fill_manual(values = c("male" = "#00BFC4", "female" = "#F8766D")) +  # swap colors
  scale_color_manual(values = c("male" = "#00BFC4", "female" = "#F8766D")) +  # swap colors
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.8),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15, colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "right"  # Ensure legend is displayed
  )
Sub.info
# figure name
fig.dpi <- 300; fig.wid <- 13; fig.hei <- 10;fig.fmt <- "png" 
fig.name <- paste("plot_subj_age_task_oldest_sex_test", ".", fig.fmt, sep = "")
## save figure
ggsave(fig.name, path = fig.savedir, Sub.info, width = fig.wid, height = fig.hei, 
       units = "cm", dpi = fig.dpi, bg = "white")

