suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/external_bulk/processed")
figure_dir <- file.path(project_root, "results/figures_external_bulk")
table_dir <- file.path(project_root, "results/tables_external_bulk")
tmp_dir <- file.path(project_root, "tmp")
for (d in c(figure_dir, table_dir, tmp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

expr <- readRDS(file.path(processed_dir, "GSE147383_expr_gene_symbol_collapsed.rds"))
pheno <- readr::read_csv(file.path(table_dir, "GSE147383_pheno_annotated.csv"), show_col_types = FALSE)
np_pheno <- pheno %>% filter(tissue == "Nucleus pulposus", !is.na(degeneration_group)) %>%
  mutate(degeneration_group = factor(degeneration_group, levels = c("young_non_degenerated", "aged_degenerated")))

selected_genes <- c("ACAN", "COL1A1", "COL3A1", "SPARC", "GPX4", "GPX3", "FTL", "SLC7A11")
selected_genes <- intersect(selected_genes, rownames(expr))

long_gene <- expr[selected_genes, np_pheno$geo_accession, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "geo_accession", values_to = "expression") %>%
  left_join(np_pheno %>% select(geo_accession, degeneration_group, age = characteristics_1_age, gender = characteristics_2_gender), by = "geo_accession")

stats <- long_gene %>%
  group_by(gene) %>%
  summarise(
    mean_young = mean(expression[degeneration_group == "young_non_degenerated"]),
    mean_aged = mean(expression[degeneration_group == "aged_degenerated"]),
    delta_aged_minus_young = mean_aged - mean_young,
    p_value = wilcox.test(expression[degeneration_group == "aged_degenerated"], expression[degeneration_group == "young_non_degenerated"], exact = FALSE)$p.value,
    .groups = "drop"
  ) %>% arrange(desc(abs(delta_aged_minus_young)))
readr::write_csv(stats, file.path(table_dir, "GSE147383_NP_supportive_selected_gene_stats.csv"))

p_gene <- ggplot(long_gene, aes(x = degeneration_group, y = expression, fill = degeneration_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.70, linewidth = 0.3) +
  geom_jitter(width = 0.08, size = 2.2, alpha = 0.9) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(young_non_degenerated = "#4E79A7", aged_degenerated = "#E15759")) +
  labs(
    title = "Supportive trend check in GSE147383 NP samples",
    subtitle = "Only 2 young non-degenerated and 2 aged degenerated NP samples; interpret direction only",
    x = NULL,
    y = "Expression (series matrix value)",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
ggsave(file.path(figure_dir, "33_GSE147383_NP_supportive_selected_gene_boxplots.png"), p_gene, width = 11.5, height = 6.4, dpi = 300)

module_sets <- list(
  Matrix_homeostasis = c("ACAN", "COL2A1", "COL6A2", "BGN"),
  Iron_ferroptosis_response = c("GPX4", "SLC39A14", "FTH1", "FTL", "SLC7A11"),
  Fibrotic_remodeling = c("COL1A1", "COL3A1", "POSTN", "ADAMTS5", "MMP13", "SPARC")
)

calc_module_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

module_df <- map_dfr(names(module_sets), function(module_name) {
  tibble(
    geo_accession = np_pheno$geo_accession,
    module = module_name,
    score = calc_module_score(expr[, np_pheno$geo_accession, drop = FALSE], module_sets[[module_name]]),
    genes_used = paste(intersect(module_sets[[module_name]], rownames(expr)), collapse = ";")
  )
}) %>% left_join(np_pheno %>% select(geo_accession, degeneration_group), by = "geo_accession")
readr::write_csv(module_df, file.path(table_dir, "GSE147383_NP_supportive_module_scores.csv"))

module_stats <- module_df %>%
  group_by(module) %>%
  summarise(
    genes_used = first(genes_used),
    mean_young = mean(score[degeneration_group == "young_non_degenerated"]),
    mean_aged = mean(score[degeneration_group == "aged_degenerated"]),
    delta_aged_minus_young = mean_aged - mean_young,
    p_value = wilcox.test(score[degeneration_group == "aged_degenerated"], score[degeneration_group == "young_non_degenerated"], exact = FALSE)$p.value,
    .groups = "drop"
  )
readr::write_csv(module_stats, file.path(table_dir, "GSE147383_NP_supportive_module_score_stats.csv"))

p_module <- ggplot(module_df, aes(x = degeneration_group, y = score, fill = degeneration_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.70, linewidth = 0.3) +
  geom_jitter(width = 0.08, size = 2.3, alpha = 0.9) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(young_non_degenerated = "#4E79A7", aged_degenerated = "#E15759")) +
  labs(
    title = "Supportive module-score trend check in GSE147383 NP samples",
    subtitle = "Small sample size: n=2 vs n=2; use only as direction-level support",
    x = NULL,
    y = "Module score",
    fill = "Group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
ggsave(file.path(figure_dir, "34_GSE147383_NP_supportive_module_score_boxplots.png"), p_module, width = 10.5, height = 4.8, dpi = 300)

print(stats, n = Inf)
print(module_stats, n = Inf)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("GSE147383 supportive plots finished.")
