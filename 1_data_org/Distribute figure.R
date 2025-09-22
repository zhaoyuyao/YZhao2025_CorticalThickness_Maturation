#------------------------------
# Project: EBDS
# Title: Subject distribution figure for EBDS data
# Author: Yuyao Zhao
# Date: Jul 17 2025
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
dir.create(fig.savedir)


# Loading packages 
packages =  c("ggplot2",'ggpubr','dplyr')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), repos = "http://cran.us.r-project.org")
}
invisible(lapply(packages, library, character.only = TRUE))

# Loading files
combined_df <- read.csv("CT_only_CS_T1_no_MVM_AAL_globZ2_train.csv")

# Reorder subjects based on first visit year, then by last visit year, and assign unique Age_ID
subject_id_df <- combined_df %>%
  group_by(Case) %>%
  summarize(First_Visit_Year = min(Year),  # Find the first visit year
            Last_Visit_Year = max(Year)) %>%  # Find the last visit year
  ungroup() %>%
  arrange(First_Visit_Year, Last_Visit_Year, Case) %>%  # Sort by first visit year, last visit year, and Subject ID
  mutate(Age_ID = row_number())  # Assign unique sequential Age_ID per subject

# Merge the Age_ID back into the original dataframe
df <- combined_df %>%
  left_join(subject_id_df %>% select(Case, Age_ID), by = "Case")

#==========Subject information============
Years <- factor(combined_df$Year) 
Timepoints <- factor(combined_df$Year)

colfunc <- colorRampPalette(c("dodgerblue4", "dodgerblue"))
colors <- colfunc(7)
plot.age <- ggplot(data = df, aes(x = Years, y = Age_ID))+
  geom_line(aes(group=Age_ID),color="gray", size=0.05)+
  geom_point(aes(group = Years, colour=Years)) + scale_color_manual(values = colors)+ 
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
  labs(x="Age (Years)",y="#Subject")
plot.age


# figure information set up
fig.dpi <- 300; fig.wid <- 13; fig.hei <- 10; fig.fmt <- "png" 
# figure name
fig.name <- paste("plot_subj_age_AAL_0to10_train", ".", fig.fmt, sep = "")
## save figure
ggsave(fig.name, path = fig.savedir, plot.age, width = fig.wid, height = fig.hei, 
       units = "cm", dpi = fig.dpi, bg = "white")


####====Sex=========
# Reorder subjects based on first visit year, then by last visit year, and assign unique Age_ID
subject_id_df <- combined_df %>%
  group_by(Case) %>%
  summarize(First_Visit_Year = min(Year),  # Find the first visit year
            Last_Visit_Year = max(Year),
            Sex = min(Sex)) %>%  # Find the last visit year
  ungroup() %>%
  arrange(Sex, First_Visit_Year, Last_Visit_Year, Case) %>%  # Sort by first visit year, last visit year, and Subject ID
  mutate(Age_ID = row_number())  # Assign unique sequential Age_ID per subject

# Merge the Age_ID back into the original dataframe
df.g <- combined_df %>%
  left_join(subject_id_df %>% select(Case, Age_ID), by = "Case")
# Change as a factor
df.g$Sex<-factor(df.g$Sex,
               levels = c('0','1'),
               labels = c("female","male"))
# plot
Sub.info <- ggplot(data = df.g, aes(x = Year, y = Age_ID, color = Sex)) +
  geom_line(aes(group = Age_ID), size = 0.05) +  # Ensure only group is set here
  geom_point() + 
  scale_x_continuous(limits = c(-0.5, 10.5), 
                     expand = c(0,0), 
                     breaks = c(0,1,2,4,6,8,10)) +
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
  labs(x="Age (Years)",y="Subject")
Sub.info

# figure name
fig.name <- paste0("plot_subj_age_AAL_0to10_sex_train", ".", fig.fmt)
fig.dpi <- 300; fig.wid <- 13; fig.hei <- 10; fig.fmt <- "png" 
## save figure
ggsave(fig.name, path = fig.savedir, Sub.info, width = fig.wid, height = fig.hei, 
       units = "cm", dpi = fig.dpi, bg = "white")
