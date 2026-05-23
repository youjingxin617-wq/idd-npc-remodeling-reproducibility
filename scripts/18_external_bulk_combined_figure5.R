suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

project_root <- normalizePath(".", winslash = "/")
figure_dir <- file.path(project_root, "results/figures_external_bulk")
table_dir <- file.path(project_root, "results/tables_external_bulk")
tmp_dir <- file.path(project_root, "tmp")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

pal <- c(
  low = "#4E79A7",
  high = "#E15759",
  young_non_degenerated = "#4E79A7",
  aged_degenerated = "#E15759"
)

module_order <- c("Matrix_homeostasis", "Iron_ferroptosis_response", "Fibrotic_remodeling")
module_labels <- c(
  Matrix_homeostasis = "Matrix\nhomeostasis",
  Iron_ferroptosis_response = "Iron/ferroptosis\nresponse",
  Fibrotic_remodeling = "Fibrotic\nremodeling"
)

gse70362_module <- readr::read_csv(file.path(table_dir, "GSE70362_NP_low_high_module_scores.csv"), show_col_types = FALSE) %>%
  mutate(
    dataset = "GSE70362",
    module = factor(module, levels = module_order),
    group = factor(group, levels = c("low", "high"))
  )
gse70362_stats <- readr::read_csv(file.path(table_dir, "GSE70362_NP_low_high_module_score_stats.csv"), show_col_types = FALSE) %>%
  mutate(module = factor(module, levels = module_order))
gse70362_auc <- readr::read_csv(file.path(table_dir, "GSE70362_NP_low_high_module_score_auc.csv"), show_col_types = FALSE) %>%
  mutate(module = factor(module, levels = module_order))

gse147383_module <- readr::read_csv(file.path(table_dir, "GSE147383_NP_supportive_module_scores.csv"), show_col_types = FALSE) %>%
  mutate(
    dataset = "GSE147383",
    module = factor(module, levels = module_order),
    degeneration_group = factor(degeneration_group, levels = c("young_non_degenerated", "aged_degenerated"))
  )
gse147383_stats <- readr::read_csv(file.path(table_dir, "GSE147383_NP_supportive_module_score_stats.csv"), show_col_types = FALSE) %>%
  mutate(module = factor(module, levels = module_order))

label_70362 <- gse70362_module %>%
  group_by(module) %>%
  summarise(y = max(score, na.rm = TRUE) + 0.12 * diff(range(score, na.rm = TRUE)), .groups = "drop") %>%
  left_join(gse70362_stats %>% select(module, p_value, p_adj), by = "module") %>%
  mutate(label = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    p_value < 0.05 ~ paste0("p=", signif(p_value, 2)),
    TRUE ~ "ns"
  ))

p_a <- ggplot(gse70362_module, aes(x = group, y = score, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.78, linewidth = 0.3) +
  geom_jitter(width = 0.10, size = 1.8, alpha = 0.85) +
  geom_text(data = label_70362, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.2) +
  facet_wrap(~ module, scales = "free_y", nrow = 1, labeller = labeller(module = module_labels)) +
  scale_fill_manual(values = pal[c("low", "high")], labels = c(low = "Low", high = "High")) +
  labs(
    title = "A. GSE70362 NP module scores",
    subtitle = "Low: Thompson I/I-II/II (n=8); High: IV/V (n=10); Grade III excluded",
    x = NULL,
    y = "Module score",
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "top"
  )

# Reconstruct ROC curves from module score table.
roc_curve <- function(score, group, label) {
  ok <- !is.na(score) & group %in% c("low", "high")
  score <- score[ok]
  group <- group[ok]
  thresholds <- sort(unique(c(Inf, score, -Inf)), decreasing = TRUE)
  map_dfr(thresholds, function(th) {
    pred <- score >= th
    pos <- group == "high"
    tibble(
      module = label,
      threshold = th,
      fpr = sum(pred & !pos) / sum(!pos),
      tpr = sum(pred & pos) / sum(pos)
    )
  })
}
roc_df <- map_dfr(unique(as.character(gse70362_module$module)), function(m) {
  dat <- gse70362_module %>% filter(as.character(module) == m)
  roc_curve(dat$score, dat$group, m)
}) %>%
  left_join(gse70362_auc, by = "module") %>%
  mutate(
    module = factor(module, levels = module_order),
    label = paste0(module_labels[as.character(module)] %>% str_replace_all("\\n", " "), " (AUC=", sprintf("%.2f", auc_direction_adjusted), ")")
  )

p_b <- ggplot(roc_df, aes(x = fpr, y = tpr, color = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.35) +
  geom_line(linewidth = 0.9) +
  coord_equal() +
  labs(
    title = "B. GSE70362 exploratory ROC",
    subtitle = "Direction-adjusted AUC; not a diagnostic model",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 7)
  )

p_c <- ggplot(gse147383_module, aes(x = degeneration_group, y = score, fill = degeneration_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.78, linewidth = 0.3) +
  geom_jitter(width = 0.07, size = 1.9, alpha = 0.9) +
  facet_wrap(~ module, scales = "free_y", nrow = 1, labeller = labeller(module = module_labels)) +
  scale_fill_manual(values = pal[c("young_non_degenerated", "aged_degenerated")], labels = c(young_non_degenerated = "Young/non-deg", aged_degenerated = "Aged/deg")) +
  labs(
    title = "C. GSE147383 supportive NP trend",
    subtitle = "Small cohort: n=2 vs n=2; direction-level support only",
    x = NULL,
    y = "Module score",
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 7),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "top",
    legend.text = element_text(size = 7)
  )

p_c_no_roc <- p_c +
  labs(title = "B. GSE147383 supportive NP trend")

combined <- (p_a / (p_b | p_c)) +
  plot_annotation(
    title = "External bulk transcriptomic support for matrix and ferroptosis-response modules",
    subtitle = "Exploratory support from independent NP transcriptomic datasets; conclusions should remain hypothesis-generating",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

ggsave(file.path(figure_dir, "35_external_bulk_combined_module_validation_Figure5.png"), combined, width = 13.5, height = 10, dpi = 300)

# A no-ROC version may be cleaner for the main manuscript.
combined_no_roc <- (p_a / p_c_no_roc) +
  plot_annotation(
    title = "External bulk transcriptomic support for matrix and ferroptosis-response module scores",
    subtitle = "GSE70362 provides exploratory support; GSE147383 provides small-cohort directional support",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )
ggsave(file.path(figure_dir, "36_external_bulk_combined_module_validation_no_ROC.png"), combined_no_roc, width = 13.5, height = 7.2, dpi = 300)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("Combined external bulk Figure 5 candidates finished.")
