# ---------------------------------------------
# CPM Observed vs Predicted (LOOCV + Test)
# ---------------------------------------------
# Usage:
#   plot_cpm("composite_score", 0.050,
#            pred_dir = "/path/to/.../predictions",
#            out_dir  = "/path/to/save/figs")
# ---------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(scales)
})

# Pretty p-value for annotation
.pp_lab <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 0.0001")          # no scientific
  paste0("p = ", formatC(p, format = "f", digits = 3))  # e.g., 0.012
}


# Core plotting function
plot_cpm <- function(behavior,
                     threshold,
                     pred_dir,
                     out_dir = pred_dir,
                     xlab = "Observed value",
                     ylab = "Predicted value",
                     xlim = NULL, ylim = NULL,
                     point_size = 2.5,
                     line_size = 0.9,
                     width = 8, height = 5.5, dpi = 300) {
  
  thr_tag   <- sprintf("thr%.3f", threshold)
  file_train <- file.path(pred_dir, sprintf("train_LOOCV_%s_%s.csv", behavior, thr_tag))
  file_test  <- file.path(pred_dir, sprintf("test_FULL_%s_%s.csv",  behavior, thr_tag))
  if (!file.exists(file_train)) stop("Missing file: ", file_train)
  if (!file.exists(file_test))  stop("Missing file: ", file_test)
  
  tr <- readr::read_csv(file_train, show_col_types = FALSE)
  te <- readr::read_csv(file_test,  show_col_types = FALSE)
  
  longify <- function(df) {
    df %>%
      mutate(Set = as.character(Set)) %>%
      tidyr::pivot_longer(
        cols = c(Pred_Combined, Pred_Pos, Pred_Neg),
        names_to = "Model", values_to = "Pred"
      ) %>%
      mutate(
        Model = dplyr::recode(Model,
                              Pred_Pos = "Positive",
                              Pred_Neg = "Negative",
                              Pred_Combined = "GLM"),
        Set = dplyr::recode(Set, `Train-LOOCV` = "Train (LOOCV)", Test = "Test")
      ) %>%
      dplyr::select(SubjectID, Set, Model, Observed, Pred, Sex, Year)
  }
  
  dat <- dplyr::bind_rows(longify(tr), longify(te)) %>%
    dplyr::filter(!is.na(Observed), !is.na(Pred)) %>%
    dplyr::mutate(
      Model = factor(Model, levels = c("Positive", "Negative", "GLM")),
      Set   = factor(Set,   levels = c("Train (LOOCV)", "Test"))
    )
  
  # Default limits if none provided (rounded to integers)
  if (is.null(xlim)) {
    xlim <- c(floor(min(dat$Observed, na.rm = TRUE)),
              ceiling(max(dat$Observed, na.rm = TRUE)))
  }
  if (is.null(ylim)) {
    ylim <- c(floor(min(dat$Pred, na.rm = TRUE)),
              ceiling(max(dat$Pred, na.rm = TRUE)))
  }
  
  # Integer tick marks based on limits
  x_breaks <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 2)
  y_breaks <- seq(floor(ylim[1]), ceiling(ylim[2]), by = 2)
  
  cols <- c("Positive" = "firebrick3", "Negative" = "dodgerblue3", "GLM" = "grey40")
  
  .pp_lab <- function(p) {
    if (is.na(p)) return("P = NA")
    if (p < 1e-4) "P < 1e-4" else paste0("p = ", format(p, digits = 2, scientific = TRUE))
  }
  
  ann <- dat %>%
    dplyr::group_by(Set, Model) %>%
    dplyr::summarise(
      r = suppressWarnings(cor(Observed, Pred, method = "spearman", use = "complete.obs")),
      p = tryCatch(suppressWarnings(cor.test(Observed, Pred, method = "spearman"))$p.value,
                   error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # italic(r) == 0.42   or   italic(r)=="NA"
      lab1 = ifelse(is.na(r), 'italic(r)=="NA"', sprintf('italic(r)==%.3f', r)),
      # italic(p) < 0.0001  or   italic(p) == 0.012  or   italic(p)=="NA"
      lab2 = dplyr::case_when(
        is.na(p)        ~ 'italic(p)=="NA"',
        p < 1e-4        ~ 'italic(p)<0.0001',
        TRUE            ~ paste0('italic(p)==', formatC(p, format = "f", digits = 3))
      )
    )
  
  x_pad <- function(v) diff(range(v, na.rm = TRUE)) * 0.02
  y_pad <- function(v) diff(range(v, na.rm = TRUE)) * 0.06
  
  pos_df <- dat %>%
    dplyr::group_by(Set, Model) %>%
    dplyr::summarise(xmin = xlim[1] + x_pad(Observed),
                     ymax = ylim[2] - y_pad(Pred),
                     .groups = "drop")
  
  ann_pos <- dplyr::left_join(ann, pos_df, by = c("Set","Model"))
  
  p <- ggplot(dat, aes(Observed, Pred, color = Model)) +
    geom_point(size = point_size, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, size = line_size, formula = y ~ x) +
    facet_grid(Set ~ Model, scales = "free_y")+
    scale_color_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = xlab, y = ylab,
         title = paste0(behavior, "  (", thr_tag, ")")) +
    theme_classic(base_size = 12) +   # axes with lines & ticks; no grid
    theme(
      # facet labels: no boxes, just text
      strip.background = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      # keep it tidy
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      # make sure ticks are visible on all panels
      axis.ticks.length = unit(0.12, "cm")
    )+
    geom_text(data = ann_pos, aes(x = xmin, y = ymax, label = lab1),
              inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 3.6, parse = TRUE) +
    geom_text(data = ann_pos, aes(x = xmin, y = ymax - y_pad(dat$Pred), label = lab2),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.6, parse = TRUE)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  base <- file.path(out_dir, sprintf("CPM_scatter_%s_%s", behavior, thr_tag))
  ggsave(paste0(base, ".pdf"), p, width = width, height = height, dpi = dpi)
  ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = dpi)
  p
}


# -------------------------
# Calls:
# -------------------------
pred_dir <- ".../3_sbMCN/sbMCN_CT_CPM/Outputs_Y0to2/reg_gender_ROIave_no_MVM_WB_AAL_Y0to2_CPM_092425/predictions"
fig_dir <- ".../3_sbMCN/sbMCN_CT_CPM/R_fig" 
plot_cpm("composite_score", 0.050, pred_dir, fig_dir)

