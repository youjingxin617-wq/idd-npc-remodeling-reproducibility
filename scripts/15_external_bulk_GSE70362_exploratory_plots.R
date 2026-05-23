suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/external_bulk/processed")
figure_dir <- file.path(project_root, "results/figures_external_bulk")
table_dir <- file.path(project_root, "results/tables_external_bulk")
tmp_dir <- file.path(project_root, "tmp")

for (d in c(figure_dir, table_dir, tmp_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

expr <- readRDS(file.path(processed_dir, "GSE70362_expr_gene_symbol_collapsed.rds"))
pheno <- readr::read_csv(file.path(table_dir, "GSE70362_pheno_annotated.csv"), show_col_types = FALSE) %>%
  mutate(
    tissue = characteristics_3_tissue,
    thompson_grade = characteristics_2_thompson.grade,
    group = case_when(
      tissue == "Nucleus pulposus" & thompson_grade %in% c("I", "I-II", "II") ~ "low",
      tissue == "Nucleus pulposus" & thompson_grade %in% c("IV", "V") ~ "high",
      tissue == "Nucleus pulposus" & thompson_grade == "III" ~ "intermediate",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("low", "high", "intermediate"))
  )

usable <- pheno %>% filter(tissue == "Nucleus pulposus", group %in% c("low", "high"))
expr_use <- expr[, usable$geo_accession, drop = FALSE]

selected_genes <- c(
  "ACAN", "COL2A1", "COL6A2", "BGN",
  "GPX4", "SLC39A14", "FTH1", "FTL",
  "COL1A1", "POSTN", "ADAMTS5", "GPX3"
)
selected_genes <- intersect(selected_genes, rownames(expr_use))

long_gene <- expr_use[selected_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "geo_accession", values_to = "expression") %>%
  left_join(usable %>% select(geo_accession, group, thompson_grade), by = "geo_accession")

stat_gene <- long_gene %>%
  group_by(gene) %>%
  summarise(
    n_low = sum(group == "low"),
    n_high = sum(group == "high"),
    mean_low = mean(expression[group == "low"]),
    mean_high = mean(expression[group == "high"]),
    delta_high_minus_low = mean_high - mean_low,
    p_value = wilcox.test(expression[group == "high"], expression[group == "low"], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_value < 0.05 ~ paste0("p=", signif(p_value, 2), " (trend)"),
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_value)
readr::write_csv(stat_gene, file.path(table_dir, "GSE70362_NP_low_high_selected_gene_plot_stats.csv"))

label_gene <- stat_gene %>%
  left_join(long_gene %>% group_by(gene) %>% summarise(y = max(expression) + 0.08 * diff(range(expression)), .groups = "drop"), by = "gene")

p_gene <- ggplot(long_gene, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.8) +
  geom_text(data = label_gene, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(low = "#4E79A7", high = "#E15759")) +
  labs(
    title = "Exploratory external validation in GSE70362 NP samples",
    subtitle = "Low degeneration: Thompson I/I-II/II; high degeneration: Thompson IV/V; Grade III excluded",
    x = NULL,
    y = "Expression (series matrix value)",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
ggsave(file.path(figure_dir, "30_GSE70362_NP_low_high_selected_gene_boxplots.png"), p_gene, width = 12, height = 8, dpi = 300)

module_sets <- list(
  Matrix_homeostasis = c("ACAN", "COL2A1", "COL6A2", "BGN"),
  Iron_ferroptosis_response = c("GPX4", "SLC39A14", "FTH1", "FTL", "SLC7A11"),
  Fibrotic_remodeling = c("COL1A1", "COL3A1", "POSTN", "ADAMTS5", "MMP13")
)

calc_module_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) < 2) {
    return(rep(NA_real_, ncol(mat)))
  }
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

module_df <- map_dfr(names(module_sets), function(module_name) {
  tibble(
    geo_accession = colnames(expr_use),
    module = module_name,
    score = calc_module_score(expr_use, module_sets[[module_name]]),
    genes_used = paste(intersect(module_sets[[module_name]], rownames(expr_use)), collapse = ";")
  )
}) %>%
  left_join(usable %>% select(geo_accession, group, thompson_grade), by = "geo_accession")
readr::write_csv(module_df, file.path(table_dir, "GSE70362_NP_low_high_module_scores.csv"))

stat_module <- module_df %>%
  filter(!is.na(score)) %>%
  group_by(module) %>%
  summarise(
    genes_used = first(genes_used),
    n_low = sum(group == "low"),
    n_high = sum(group == "high"),
    mean_low = mean(score[group == "low"]),
    mean_high = mean(score[group == "high"]),
    delta_high_minus_low = mean_high - mean_low,
    p_value = wilcox.test(score[group == "high"], score[group == "low"], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_value < 0.05 ~ paste0("p=", signif(p_value, 2), " (trend)"),
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_value)
readr::write_csv(stat_module, file.path(table_dir, "GSE70362_NP_low_high_module_score_stats.csv"))

label_module <- stat_module %>%
  left_join(module_df %>% group_by(module) %>% summarise(y = max(score, na.rm = TRUE) + 0.10 * diff(range(score, na.rm = TRUE)), .groups = "drop"), by = "module")

p_module <- ggplot(module_df, aes(x = group, y = score, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 2.0, alpha = 0.85) +
  geom_text(data = label_module, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(low = "#4E79A7", high = "#E15759")) +
  labs(
    title = "Exploratory module-score validation in GSE70362 NP samples",
    subtitle = "Scores are averaged z-scores of genes available in the expression matrix",
    x = NULL,
    y = "Module score",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
ggsave(file.path(figure_dir, "31_GSE70362_NP_low_high_module_score_boxplots.png"), p_module, width = 10.5, height = 4.8, dpi = 300)

# Simple ROC/AUC calculation without external packages. Positive class = high.
auc_rank <- function(score, group) {
  ok <- !is.na(score) & group %in% c("low", "high")
  score <- score[ok]
  group <- group[ok]
  y <- as.integer(group == "high")
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  auc <- (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  auc
}

roc_curve <- function(score, group, label) {
  ok <- !is.na(score) & group %in% c("low", "high")
  score <- score[ok]
  group <- group[ok]
  thresholds <- sort(unique(c(Inf, score, -Inf)), decreasing = TRUE)
  map_dfr(thresholds, function(th) {
    pred <- score >= th
    pos <- group == "high"
    tibble(
      label = label,
      threshold = th,
      fpr = sum(pred & !pos) / sum(!pos),
      tpr = sum(pred & pos) / sum(pos)
    )
  })
}

roc_modules <- module_df %>%
  group_by(module) %>%
  summarise(
    auc_high_when_score_higher = auc_rank(score, group),
    auc_direction_adjusted = pmax(auc_high_when_score_higher, 1 - auc_high_when_score_higher),
    direction = if_else(auc_high_when_score_higher >= 0.5, "higher_score_in_high", "lower_score_in_high"),
    .groups = "drop"
  )
readr::write_csv(roc_modules, file.path(table_dir, "GSE70362_NP_low_high_module_score_auc.csv"))

roc_df <- map_dfr(unique(module_df$module), function(m) {
  dat <- module_df %>% filter(module == m)
  roc_curve(dat$score, dat$group, m)
}) %>%
  left_join(roc_modules, by = c("label" = "module")) %>%
  mutate(label_auc = paste0(label, " (AUC=", sprintf("%.2f", auc_direction_adjusted), ")"))

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = label_auc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888") +
  geom_line(linewidth = 0.9) +
  coord_equal() +
  labs(
    title = "Exploratory ROC of module scores in GSE70362 NP samples",
    subtitle = "AUC is direction-adjusted because some biological scores decrease in high degeneration",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "bottom"
  )
ggsave(file.path(figure_dir, "32_GSE70362_NP_low_high_module_score_ROC.png"), p_roc, width = 7.2, height = 6.2, dpi = 300)

print(stat_gene, n = Inf)
print(stat_module, n = Inf)
print(roc_modules, n = Inf)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("GSE70362 exploratory external validation plots finished.")
