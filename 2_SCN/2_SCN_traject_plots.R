#------------------------------
# Project: EBDS
# Title: SCN trajectories graph
# Author: Yuyao Zhao
# Date: Aug 5 2025
#------------------------------
rm(list = ls())

library(ggplot2)
library(dplyr)

# --- Paths ---
setwd('...')
cur.dir <- getwd()
fig.fold <- "R_figures"                       # fix: spelling
fig.savedir <- file.path(cur.dir, fig.fold)
dir.create(fig.savedir, showWarnings = FALSE, recursive = TRUE)

# --- Data ---
df <- read.csv("SCN_metrics_long_globZ2_train.csv")
# ensure expected columns exist
stopifnot(all(c("Measure","Year","Network","Mean","STD","p_value") %in% names(df)))
df$Year.f <- factor(df$Year)

# p -> stars
p_to_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
df <- df %>% mutate(sig_label = p_to_stars(p_value))

# Colors
network_colors <- c(
  "Vis" = 'dodgerblue',
  "Som" = 'chartreuse4',
  "Lim" = 'coral',
  "Pos" = 'firebrick3',
  "Neg" = 'dodgerblue3',
  "WB"  = 'gray'
)

# Year-pair map for bottom annotations (0-1, 1-2, 2-4, 4-6, 6-8, 8-10)
year_pairs <- data.frame(
  Year     = c(0, 1, 2, 4, 6, 8),
  Year_end = c(1, 2, 4, 6, 8, 10)
)

# ---- Plot function ----
plot_metric <- function(metric_name, y_lab = NULL, y_min = NULL) {
  df_metric <- df %>% filter(Measure == metric_name)
  
  # Order networks by mean at Year==0 (bottom to top)
  network_order <- df_metric %>%
    filter(Year == 0) %>%
    arrange(Mean) %>%
    pull(Network)
  network_order <- unique(network_order)  # safety
  df_metric$Network <- factor(df_metric$Network, levels = network_order)
  
  # Bottom annotation frame
  df_annot <- df_metric %>%
    left_join(year_pairs, by = "Year") %>%
    filter(!is.na(sig_label), sig_label != "", !is.na(Year_end)) %>%
    mutate(
      net_index = as.integer(factor(Network, levels = network_order))
    )
  
  # Compute a consistent bottom band so lines + stars don’t overlap the data
  # Use data range incl. error bars to place the bottom band a bit lower.
  y_min_data <- min(df_metric$Mean - df_metric$STD, na.rm = TRUE)
  y_max_data <- max(df_metric$Mean + df_metric$STD, na.rm = TRUE)
  y_range    <- y_max_data - y_min_data
  band_floor <- y_min_data - 0.22 * y_range              # a bit below data
  step       <- 0.04 * y_range                            # spacing between nets
  
  df_annot <- df_annot %>%
    mutate(y_pos = band_floor + (net_index - 1) * step)
  
  # X scale (unique sorted breaks)
  x_breaks <- sort(unique(df_metric$Year))
  
  # Y label default
  if (is.null(y_lab)) y_lab <- metric_name
  
  p <- ggplot(df_metric, aes(x = Year, y = Mean, color = Network, group = Network)) +
    # error bars using STD
    geom_errorbar(aes(ymin = Mean - STD, ymax = Mean + STD), width = 0.2, alpha = 0.8) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = x_breaks, labels = x_breaks,
                       expand = expansion(mult = c(0.02, 0.12))) +  # extra bottom room for annotations
    scale_color_manual(values = network_colors) +
    labs(x = "Age (Years)", y = y_lab) +
    theme(
      panel.grid.major = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.8),
      axis.line.y = element_line(colour = "black", linewidth = 0.8),
      axis.text = element_text(size = 13, colour = "black"),
      axis.title = element_text(size = 15, colour = "black"),
      panel.background = element_rect(fill = "transparent"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent"),
      legend.position = ""
    )
  
  # Optional fixed lower limit (useful if you want Degree to start at 0)
  if (!is.null(y_min)) {
    p <- p + coord_cartesian(ylim = c(y_min, NA), clip = "on")
  }
  
  # Bottom significance lines + stars
  p <- p +
    geom_segment(
      data = df_annot,
      aes(x = Year, xend = Year_end, y = y_pos, yend = y_pos, color = Network),
      inherit.aes = FALSE,
      linewidth = 0.5
    ) +
    geom_text(
      data = df_annot,
      aes(x = (Year + Year_end) / 2, y = y_pos, label = sig_label, color = Network),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0.3
    )
  
  p
}

# ---- Build plots ----
p_Seg  <- plot_metric("Seg",  y_lab = "System Segregation")
p_GE   <- plot_metric("GE",   y_lab = "Global Efficiency")
p_Deg  <- plot_metric("Degree", y_lab = "Average Degree", y_min = 0)  # start from 0 if you like

# Print to viewer
p_Seg; p_GE; p_MOD; p_Deg; p_Betw

# ---- Save figures ----
fig.wid <- 11; fig.hei <- 10; fig.dpi <- 300

ggsave(filename = "Segregation_plot.png",       plot = p_Seg,  path = fig.savedir,
       width = fig.wid, height = fig.hei, units = "cm", dpi = fig.dpi, bg = "white")
ggsave(filename = "Global_Efficiency_plot.png", plot = p_GE,   path = fig.savedir,
       width = fig.wid, height = fig.hei, units = "cm", dpi = fig.dpi, bg = "white")
ggsave(filename = "Average_Degree_plot.png",    plot = p_Deg,  path = fig.savedir,
       width = fig.wid, height = fig.hei, units = "cm", dpi = fig.dpi, bg = "white")

