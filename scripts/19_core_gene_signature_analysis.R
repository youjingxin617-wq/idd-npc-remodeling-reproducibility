suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures_core_signature")
table_dir <- file.path(project_root, "results/tables_core_signature")
bulk_table_dir <- file.path(project_root, "results/tables_external_bulk")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260521)

theme_sig <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "top"
    )
}

pal_group <- c(
  mild = "#4E79A7",
  severe = "#E15759",
  low = "#4E79A7",
  high = "#E15759",
  young_non_degenerated = "#4E79A7",
  aged_degenerated = "#E15759"
)

core_gene_sets <- list(
  Matrix_homeostasis = c("ACAN", "COL2A1", "COL6A2", "BGN"),
  Iron_ferroptosis_response = c("GPX4", "SLC39A14", "FTH1", "FTL", "SLC7A11"),
  Fibrotic_remodeling = c("COL1A1", "COL3A1", "POSTN", "ADAMTS5", "MMP13")
)

positive_genes <- c(core_gene_sets$Iron_ferroptosis_response, core_gene_sets$Fibrotic_remodeling)
negative_genes <- core_gene_sets$Matrix_homeostasis
core_genes <- unique(c(negative_genes, positive_genes))

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"

genes_found <- intersect(core_genes, rownames(clean_npc))
if (length(genes_found) < 8) {
  stop("Too few core signature genes were found in the NPC expression matrix: ", length(genes_found))
}

gene_set_table <- tibble(
  module = names(core_gene_sets),
  role_in_signature = c("negative arm", "positive arm", "positive arm"),
  genes_defined = map_chr(core_gene_sets, paste, collapse = ";"),
  genes_found = map_chr(core_gene_sets, ~ paste(intersect(.x, rownames(clean_npc)), collapse = ";")),
  n_defined = lengths(core_gene_sets),
  n_found = map_int(core_gene_sets, ~ length(intersect(.x, rownames(clean_npc))))
)
write_csv(gene_set_table, file.path(table_dir, "core_gene_signature_gene_sets.csv"))

expr <- as.matrix(GetAssayData(clean_npc, assay = "RNA", layer = "data")[genes_found, , drop = FALSE])
scaled_expr <- t(scale(t(expr)))
scaled_expr[is.na(scaled_expr)] <- 0

positive_found <- intersect(positive_genes, rownames(scaled_expr))
negative_found <- intersect(negative_genes, rownames(scaled_expr))

matrix_score <- colMeans(scaled_expr[negative_found, , drop = FALSE])
remodel_iron_score <- colMeans(scaled_expr[positive_found, , drop = FALSE])
core_signature <- remodel_iron_score - matrix_score

clean_npc$Matrix_homeostasis_core <- matrix_score
clean_npc$Remodeling_iron_core <- remodel_iron_score
clean_npc$IDD_core_signature <- core_signature

score_df <- clean_npc@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  select(
    cell_barcode,
    sample_id,
    grade,
    article_group,
    npc_subtype_final,
    Matrix_homeostasis_core,
    Remodeling_iron_core,
    IDD_core_signature
  ) %>%
  mutate(
    article_group = factor(article_group, levels = c("mild", "severe")),
    npc_subtype_final = factor(
      npc_subtype_final,
      levels = c(
        "ECM/Adh-NPCs", "Hom-NPCs", "Eff-NPCs", "Ht-NPCs",
        "Fibro-NPCs", "Fibro-reg NPCs", "Reg-NPCs"
      )
    )
  )

write_csv(score_df, file.path(table_dir, "clean_npc_core_signature_scores_per_cell.csv"))

score_summary <- score_df %>%
  pivot_longer(
    c(Matrix_homeostasis_core, Remodeling_iron_core, IDD_core_signature),
    names_to = "score",
    values_to = "value"
  ) %>%
  group_by(score, npc_subtype_final, article_group) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(value, na.rm = TRUE),
    median_score = median(value, na.rm = TRUE),
    sd_score = sd(value, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(score_summary, file.path(table_dir, "clean_npc_core_signature_summary_by_subtype_group.csv"))

score_stats <- score_df %>%
  pivot_longer(
    c(Matrix_homeostasis_core, Remodeling_iron_core, IDD_core_signature),
    names_to = "score",
    values_to = "value"
  ) %>%
  group_by(score, npc_subtype_final) %>%
  summarise(
    n_mild = sum(article_group == "mild" & !is.na(value)),
    n_severe = sum(article_group == "severe" & !is.na(value)),
    mean_mild = mean(value[article_group == "mild"], na.rm = TRUE),
    mean_severe = mean(value[article_group == "severe"], na.rm = TRUE),
    delta_severe_minus_mild = mean_severe - mean_mild,
    p_value = {
      mild_values <- value[article_group == "mild" & !is.na(value)]
      severe_values <- value[article_group == "severe" & !is.na(value)]
      if (length(mild_values) >= 3 && length(severe_values) >= 3) {
        wilcox.test(mild_values, severe_values, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  group_by(score) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()
write_csv(score_stats, file.path(table_dir, "clean_npc_core_signature_wilcox_stats.csv"))

label_df <- score_stats %>%
  filter(score == "IDD_core_signature") %>%
  select(npc_subtype_final, p_label) %>%
  left_join(
    score_df %>%
      group_by(npc_subtype_final) %>%
      summarise(y = max(IDD_core_signature, na.rm = TRUE), .groups = "drop") %>%
      mutate(y = y + 0.08 * diff(range(score_df$IDD_core_signature, na.rm = TRUE))),
    by = "npc_subtype_final"
  )

p_core_violin <- ggplot(score_df, aes(x = npc_subtype_final, y = IDD_core_signature, fill = article_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_violin(
    position = position_dodge(width = 0.82),
    scale = "width",
    trim = TRUE,
    linewidth = 0.25,
    alpha = 0.86
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.82),
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.72,
    linewidth = 0.25
  ) +
  geom_text(
    data = label_df,
    aes(x = npc_subtype_final, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  scale_fill_manual(values = pal_group[c("mild", "severe")]) +
  labs(
    title = "A. Core gene signature across NPC subtypes",
    subtitle = "IDD core signature = remodeling/iron arm - matrix-homeostasis arm",
    x = NULL,
    y = "Core gene signature score",
    fill = "Degeneration group"
  ) +
  theme_sig(10)
ggsave(file.path(figure_dir, "37_clean_npc_core_gene_signature_violin.png"), p_core_violin, width = 11, height = 5.8, dpi = 300)

p_score_components <- score_df %>%
  pivot_longer(
    c(Matrix_homeostasis_core, Remodeling_iron_core, IDD_core_signature),
    names_to = "score",
    values_to = "value"
  ) %>%
  mutate(
    score = factor(
      score,
      levels = c("Matrix_homeostasis_core", "Remodeling_iron_core", "IDD_core_signature"),
      labels = c("Matrix homeostasis", "Remodeling + iron", "IDD core signature")
    )
  ) %>%
  ggplot(aes(x = article_group, y = value, fill = article_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.78, linewidth = 0.3) +
  geom_jitter(width = 0.08, size = 0.25, alpha = 0.18) +
  facet_wrap(~ score, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_group[c("mild", "severe")]) +
  labs(
    title = "Global NPC-level core signature components",
    subtitle = "Each dot is one NPC; use as exploratory cell-level evidence",
    x = NULL,
    y = "Scaled score",
    fill = "Degeneration group"
  ) +
  theme_sig(10) +
  theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(figure_dir, "38_clean_npc_core_signature_global_components.png"), p_score_components, width = 10.5, height = 4.8, dpi = 300)

avg_expr <- score_df %>%
  mutate(subtype_group = paste(npc_subtype_final, article_group, sep = "_")) %>%
  select(cell_barcode, subtype_group) %>%
  left_join(
    as.data.frame(t(scaled_expr[genes_found, , drop = FALSE])) %>%
      rownames_to_column("cell_barcode"),
    by = "cell_barcode"
  ) %>%
  pivot_longer(all_of(genes_found), names_to = "gene", values_to = "z_expr") %>%
  group_by(gene, subtype_group) %>%
  summarise(mean_z = mean(z_expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = subtype_group, values_from = mean_z) %>%
  arrange(match(gene, core_genes))

write_csv(avg_expr, file.path(table_dir, "clean_npc_core_gene_average_scaled_expression.csv"))

heat_mat <- avg_expr %>%
  column_to_rownames("gene") %>%
  as.matrix()

png(file.path(figure_dir, "39_clean_npc_core_gene_signature_heatmap.png"), width = 2400, height = 1450, res = 220)
pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("#2D5F9A", "white", "#B53636"))(100),
  fontsize = 8,
  main = "Core signature genes across NPC subtype and degeneration groups"
)
dev.off()

p_umap <- FeaturePlot(
  clean_npc,
  features = "IDD_core_signature",
  reduction = if ("umap" %in% Reductions(clean_npc)) "umap" else "tsne",
  cols = c("#D9E7F5", "#B53636"),
  pt.size = 0.25
) +
  ggtitle("IDD core signature on NPC embedding") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(figure_dir, "40_clean_npc_core_gene_signature_embedding.png"), p_umap, width = 7.2, height = 5.8, dpi = 300)

bulk_outputs <- list()

gse70362_path <- file.path(bulk_table_dir, "GSE70362_NP_low_high_module_scores.csv")
if (file.exists(gse70362_path)) {
  gse70362_core <- read_csv(gse70362_path, show_col_types = FALSE) %>%
    select(sample_id = geo_accession, group, module, score) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    mutate(
      dataset = "GSE70362",
      comparison_group = recode(group, low = "Low", high = "High"),
      comparison_order = recode(group, low = 1L, high = 2L),
      Core_gene_signature = Iron_ferroptosis_response + Fibrotic_remodeling - Matrix_homeostasis
    ) %>%
    select(dataset, sample_id, comparison_group, comparison_order, Core_gene_signature)
  bulk_outputs$GSE70362 <- gse70362_core
}

gse147383_path <- file.path(bulk_table_dir, "GSE147383_NP_supportive_module_scores.csv")
if (file.exists(gse147383_path)) {
  gse147383_core <- read_csv(gse147383_path, show_col_types = FALSE) %>%
    select(sample_id = geo_accession, degeneration_group, module, score) %>%
    pivot_wider(names_from = module, values_from = score) %>%
    mutate(
      dataset = "GSE147383",
      comparison_group = recode(
        degeneration_group,
        young_non_degenerated = "Young/non-deg",
        aged_degenerated = "Aged/deg"
      ),
      comparison_order = recode(degeneration_group, young_non_degenerated = 1L, aged_degenerated = 2L),
      Core_gene_signature = Iron_ferroptosis_response + Fibrotic_remodeling - Matrix_homeostasis
    ) %>%
    select(dataset, sample_id, comparison_group, comparison_order, Core_gene_signature)
  bulk_outputs$GSE147383 <- gse147383_core
}

bulk_core <- bind_rows(bulk_outputs)
if (nrow(bulk_core) > 0) {
  bulk_core <- bulk_core %>%
    mutate(
      dataset = factor(dataset, levels = c("GSE70362", "GSE147383")),
      comparison_group = factor(
        comparison_group,
        levels = c("Low", "High", "Young/non-deg", "Aged/deg")
      )
    )
  write_csv(bulk_core, file.path(table_dir, "external_bulk_core_gene_signature_scores.csv"))

  bulk_stats <- bulk_core %>%
    arrange(dataset, comparison_order) %>%
    group_by(dataset) %>%
    summarise(
      group_1 = first(as.character(comparison_group[comparison_order == 1])),
      group_2 = first(as.character(comparison_group[comparison_order == 2])),
      n_group_1 = sum(comparison_group == group_1, na.rm = TRUE),
      n_group_2 = sum(comparison_group == group_2, na.rm = TRUE),
      mean_group_1 = mean(Core_gene_signature[comparison_group == group_1], na.rm = TRUE),
      mean_group_2 = mean(Core_gene_signature[comparison_group == group_2], na.rm = TRUE),
      delta_group_2_minus_group_1 = mean_group_2 - mean_group_1,
      p_value = {
        v1 <- Core_gene_signature[comparison_group == group_1]
        v2 <- Core_gene_signature[comparison_group == group_2]
        if (length(v1) >= 3 && length(v2) >= 3) {
          wilcox.test(v1, v2, exact = FALSE)$p.value
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )
  write_csv(bulk_stats, file.path(table_dir, "external_bulk_core_gene_signature_stats.csv"))

  p_bulk <- ggplot(bulk_core, aes(x = comparison_group, y = Core_gene_signature, fill = comparison_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.78, linewidth = 0.3) +
    geom_jitter(width = 0.08, size = 1.8, alpha = 0.86) +
    facet_wrap(~ dataset, scales = "free", nrow = 1) +
    scale_x_discrete(drop = TRUE) +
    scale_fill_manual(
      values = c("Low" = "#4E79A7", "High" = "#E15759", "Young/non-deg" = "#4E79A7", "Aged/deg" = "#E15759"),
      drop = TRUE
    ) +
    labs(
      title = "B. External bulk support for the core gene signature",
      subtitle = "Composite score from module scores; exploratory/supportive evidence only",
      x = NULL,
      y = "Core gene signature score",
      fill = NULL
    ) +
    theme_sig(10)
  ggsave(file.path(figure_dir, "41_external_bulk_core_gene_signature_validation.png"), p_bulk, width = 8.8, height = 4.8, dpi = 300)
}

combined <- (p_core_violin / (p_score_components | p_bulk)) +
  plot_annotation(
    title = "Core gene signature links NPC state shifts with external bulk trends",
    subtitle = "Signature increases when remodeling/iron-response genes rise relative to matrix-homeostasis genes",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )
ggsave(file.path(figure_dir, "42_core_gene_signature_combined_Figure6_candidate.png"), combined, width = 13.5, height = 10, dpi = 300)

main_combined <- (p_core_violin / p_bulk) +
  plot_annotation(
    title = "Core gene signature links NPC state shifts with external bulk trends",
    subtitle = "Signature increases when remodeling/iron-response genes rise relative to matrix-homeostasis genes",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )
ggsave(file.path(figure_dir, "43_core_gene_signature_main_Figure6_candidate.png"), main_combined, width = 12.5, height = 9, dpi = 300)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("Core gene signature analysis finished.")
