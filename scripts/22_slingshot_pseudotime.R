suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(scater)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures_slingshot")
table_dir <- file.path(project_root, "results/tables_slingshot")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260521)

subtype_levels <- c(
  "ECM/Adh-NPCs", "Hom-NPCs", "Eff-NPCs", "Ht-NPCs",
  "Fibro-NPCs", "Fibro-reg NPCs", "Hom/Reg-like NPCs"
)
subtype_pal <- c(
  "ECM/Adh-NPCs" = "#4E79A7",
  "Hom-NPCs" = "#59A14F",
  "Eff-NPCs" = "#F28E2B",
  "Ht-NPCs" = "#B07AA1",
  "Fibro-NPCs" = "#E15759",
  "Fibro-reg NPCs" = "#9C755F",
  "Hom/Reg-like NPCs" = "#76B7B2"
)
pal_group <- c(mild = "#4E79A7", severe = "#E15759")

safe_get_assay <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

theme_pt <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

normalize_pt <- function(mat) {
  out <- mat
  for (j in seq_len(ncol(out))) {
    x <- out[, j]
    rng <- range(x, na.rm = TRUE)
    if (all(is.finite(rng)) && diff(rng) > 0) {
      out[, j] <- (x - rng[1]) / diff(rng)
    }
  }
  out
}

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"

meta <- clean_npc[[]] %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    npc_subtype_final = factor(npc_subtype_final, levels = subtype_levels),
    article_group = factor(article_group, levels = c("mild", "severe"))
  ) %>%
  filter(!is.na(npc_subtype_final), !is.na(article_group))

clean_npc <- subset(clean_npc, cells = meta$cell_barcode)
meta <- clean_npc[[]] %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    npc_subtype_final = factor(npc_subtype_final, levels = subtype_levels),
    article_group = factor(article_group, levels = c("mild", "severe"))
  )

emb <- Embeddings(clean_npc, "umap")[meta$cell_barcode, 1:2, drop = FALSE]
colnames(emb) <- c("UMAP_1", "UMAP_2")
meta <- bind_cols(meta, as.data.frame(emb))

expr_data <- safe_get_assay(clean_npc, layer = "data")
expr_counts <- safe_get_assay(clean_npc, layer = "counts")
expr_data <- expr_data[, meta$cell_barcode, drop = FALSE]
expr_counts <- expr_counts[, meta$cell_barcode, drop = FALSE]

col_data <- as.data.frame(meta)
rownames(col_data) <- col_data$cell_barcode
sce <- SingleCellExperiment(
  assays = list(counts = expr_counts, logcounts = expr_data),
  colData = DataFrame(col_data)
)
reducedDims(sce)$UMAP <- as.matrix(emb)

root_subtype <- "Hom/Reg-like NPCs"
if (!root_subtype %in% as.character(meta$npc_subtype_final)) {
  root_subtype <- names(sort(table(meta$npc_subtype_final), decreasing = TRUE))[1]
}

sce <- slingshot(
  sce,
  clusterLabels = "npc_subtype_final",
  reducedDim = "UMAP",
  start.clus = root_subtype
)

lineages <- slingLineages(sce)
lineage_tbl <- tibble(
  lineage = paste0("Lineage", seq_along(lineages)),
  clusters = vapply(lineages, paste, collapse = " -> ", FUN.VALUE = character(1))
)
write_csv(lineage_tbl, file.path(table_dir, "npc_slingshot_lineages.csv"))

pt_raw <- slingPseudotime(sce, na = TRUE)
pt_norm <- normalize_pt(pt_raw)
composite_pt <- rowMeans(pt_norm, na.rm = TRUE)

if (all(is.na(composite_pt))) {
  pt_raw <- slingPseudotime(sce, na = FALSE)
  pt_norm <- normalize_pt(pt_raw)
  composite_pt <- rowMeans(pt_norm, na.rm = TRUE)
}
composite_pt <- (composite_pt - min(composite_pt, na.rm = TRUE)) / diff(range(composite_pt, na.rm = TRUE))

pt_df <- as_tibble(pt_raw, .name_repair = "unique")
colnames(pt_df) <- paste0("slingshot_raw_", seq_len(ncol(pt_df)))
pt_norm_df <- as_tibble(pt_norm, .name_repair = "unique")
colnames(pt_norm_df) <- paste0("slingshot_norm_", seq_len(ncol(pt_norm_df)))

meta_pt <- meta %>%
  bind_cols(pt_df, pt_norm_df) %>%
  mutate(slingshot_pseudotime = composite_pt)

mst_file <- file.path(project_root, "results/tables_pseudotime/npc_mst_exploratory_pseudotime_per_cell.csv")
if (file.exists(mst_file)) {
  mst_pt <- read_csv(mst_file, show_col_types = FALSE) %>%
    select(cell_barcode, mst_pseudotime = pseudotime)
  meta_pt <- left_join(meta_pt, mst_pt, by = "cell_barcode")
}

write_csv(meta_pt, file.path(table_dir, "npc_slingshot_pseudotime_per_cell.csv"))

summary_subtype <- meta_pt %>%
  group_by(npc_subtype_final, article_group) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(slingshot_pseudotime, na.rm = TRUE),
    median_pseudotime = median(slingshot_pseudotime, na.rm = TRUE),
    sd_pseudotime = sd(slingshot_pseudotime, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_subtype, file.path(table_dir, "npc_slingshot_pseudotime_summary_by_subtype_group.csv"))

stats_global <- meta_pt %>%
  summarise(
    n_mild = sum(article_group == "mild"),
    n_severe = sum(article_group == "severe"),
    mean_mild = mean(slingshot_pseudotime[article_group == "mild"], na.rm = TRUE),
    mean_severe = mean(slingshot_pseudotime[article_group == "severe"], na.rm = TRUE),
    delta_severe_minus_mild = mean_severe - mean_mild,
    p_value = wilcox.test(
      slingshot_pseudotime[article_group == "mild"],
      slingshot_pseudotime[article_group == "severe"],
      exact = FALSE
    )$p.value,
    mst_correlation = if ("mst_pseudotime" %in% colnames(meta_pt)) {
      cor(slingshot_pseudotime, mst_pseudotime, use = "complete.obs", method = "spearman")
    } else {
      NA_real_
    }
  )
write_csv(stats_global, file.path(table_dir, "npc_slingshot_pseudotime_global_mild_vs_severe_stats.csv"))

stats_subtype <- meta_pt %>%
  group_by(npc_subtype_final) %>%
  summarise(
    n_mild = sum(article_group == "mild"),
    n_severe = sum(article_group == "severe"),
    mean_mild = mean(slingshot_pseudotime[article_group == "mild"], na.rm = TRUE),
    mean_severe = mean(slingshot_pseudotime[article_group == "severe"], na.rm = TRUE),
    delta_severe_minus_mild = mean_severe - mean_mild,
    p_value = if (n_mild >= 3 && n_severe >= 3) {
      wilcox.test(
        slingshot_pseudotime[article_group == "mild"],
        slingshot_pseudotime[article_group == "severe"],
        exact = FALSE
      )$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
write_csv(stats_subtype, file.path(table_dir, "npc_slingshot_pseudotime_mild_vs_severe_stats_by_subtype.csv"))

curves <- slingCurves(sce)
curve_df <- purrr::imap_dfr(curves, function(curve, nm) {
  as_tibble(curve$s[, 1:2, drop = FALSE]) %>%
    set_names(c("UMAP_1", "UMAP_2")) %>%
    mutate(lineage = nm, order = row_number())
})

centroids <- meta_pt %>%
  group_by(npc_subtype_final) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    n_cells = n(),
    .groups = "drop"
  )
write_csv(centroids, file.path(table_dir, "npc_slingshot_subtype_centroids.csv"))

p_lineage <- ggplot(meta_pt, aes(UMAP_1, UMAP_2, color = npc_subtype_final)) +
  geom_point(size = 0.15, alpha = 0.45) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2, group = lineage), inherit.aes = FALSE, color = "black", linewidth = 0.85, alpha = 0.8, arrow = arrow(length = unit(0.09, "inches"))) +
  geom_point(data = centroids, aes(UMAP_1, UMAP_2), inherit.aes = FALSE, size = 2.5, shape = 21, fill = "white", color = "black", stroke = 0.7) +
  geom_text_repel(data = centroids, aes(UMAP_1, UMAP_2, label = npc_subtype_final), inherit.aes = FALSE, size = 2.8, max.overlaps = 20) +
  scale_color_manual(values = subtype_pal, drop = FALSE) +
  labs(
    title = "A. Slingshot-inferred NPC lineage structure",
    subtitle = paste0("Root cluster: ", root_subtype, "; curves are principal curves in UMAP space"),
    x = "UMAP 1",
    y = "UMAP 2",
    color = "NPC subtype"
  ) +
  theme_pt(9)

p_pt <- ggplot(meta_pt, aes(UMAP_1, UMAP_2, color = slingshot_pseudotime)) +
  geom_point(size = 0.15, alpha = 0.62) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2, group = lineage), inherit.aes = FALSE, color = "black", linewidth = 0.75, alpha = 0.75, arrow = arrow(length = unit(0.09, "inches"))) +
  scale_color_gradientn(colors = c("#2D5F9A", "#F4E8B5", "#B53636")) +
  labs(
    title = "B. Composite Slingshot pseudotime",
    subtitle = "Lineage pseudotimes were normalized and averaged per cell",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"
  ) +
  theme_pt(9)

p_group <- ggplot(meta_pt, aes(UMAP_1, UMAP_2, color = article_group)) +
  geom_point(size = 0.15, alpha = 0.52) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2, group = lineage), inherit.aes = FALSE, color = "black", linewidth = 0.65, alpha = 0.75, arrow = arrow(length = unit(0.09, "inches"))) +
  scale_color_manual(values = pal_group) +
  labs(
    title = "C. Degeneration groups along Slingshot trajectory",
    subtitle = "Descriptive distribution; not proof of real temporal progression",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Group"
  ) +
  theme_pt(9)

p_box_global <- ggplot(meta_pt, aes(article_group, slingshot_pseudotime, fill = article_group)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.82, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.18, alpha = 0.16) +
  scale_fill_manual(values = pal_group) +
  labs(
    title = "D. Global Slingshot pseudotime by group",
    subtitle = paste0("Severe - mild mean delta = ", signif(stats_global$delta_severe_minus_mild, 3)),
    x = NULL,
    y = "Slingshot pseudotime",
    fill = "Group"
  ) +
  theme_pt(10) +
  theme(legend.position = "none")

p_box_subtype <- ggplot(meta_pt, aes(npc_subtype_final, slingshot_pseudotime, fill = article_group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.82, linewidth = 0.3, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = article_group), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 0.16, alpha = 0.16, show.legend = FALSE) +
  scale_fill_manual(values = pal_group) +
  scale_color_manual(values = pal_group) +
  labs(
    title = "E. Slingshot pseudotime by NPC subtype and group",
    subtitle = "Cell-level comparison; interpret with subtype composition",
    x = NULL,
    y = "Slingshot pseudotime",
    fill = "Group"
  ) +
  theme_pt(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

core_gene_sets <- list(
  Matrix_homeostasis = c("ACAN", "COL2A1", "COL6A2", "BGN"),
  Iron_ferroptosis_response = c("GPX4", "SLC39A14", "FTH1", "FTL", "SLC7A11"),
  Fibrotic_remodeling = c("COL1A1", "COL3A1", "POSTN", "ADAMTS5", "MMP13")
)
trend_genes <- c("ACAN", "COL2A1", "COL1A1", "COL3A1", "POSTN", "MMP13", "GPX4", "SLC39A14", "FTH1", "FTL")
genes_use <- intersect(unique(c(trend_genes, unlist(core_gene_sets))), rownames(clean_npc))
expr_dense <- as.matrix(expr_data[genes_use, meta_pt$cell_barcode, drop = FALSE])
scaled_expr <- t(scale(t(expr_dense)))
scaled_expr[is.na(scaled_expr)] <- 0

module_mat <- purrr::map_dfc(names(core_gene_sets), function(module_name) {
  genes <- intersect(core_gene_sets[[module_name]], rownames(scaled_expr))
  tibble(!!module_name := colMeans(scaled_expr[genes, , drop = FALSE]))
})
module_df <- tibble(cell_barcode = colnames(scaled_expr)) %>%
  bind_cols(module_mat) %>%
  left_join(meta_pt %>% select(cell_barcode, slingshot_pseudotime, article_group, npc_subtype_final), by = "cell_barcode") %>%
  mutate(
    Composite_remodeling_iron_over_matrix = Iron_ferroptosis_response + Fibrotic_remodeling - Matrix_homeostasis,
    pseudotime_bin = cut(slingshot_pseudotime, breaks = seq(0, 1, length.out = 26), include.lowest = TRUE, labels = FALSE)
  )
write_csv(module_df, file.path(table_dir, "npc_slingshot_module_scores_per_cell.csv"))

module_trend_df <- module_df %>%
  pivot_longer(c(Matrix_homeostasis, Iron_ferroptosis_response, Fibrotic_remodeling, Composite_remodeling_iron_over_matrix), names_to = "module", values_to = "score") %>%
  group_by(module, pseudotime_bin) %>%
  summarise(
    pseudotime = mean(slingshot_pseudotime, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    module = factor(
      module,
      levels = c("Matrix_homeostasis", "Iron_ferroptosis_response", "Fibrotic_remodeling", "Composite_remodeling_iron_over_matrix"),
      labels = c("Matrix", "Iron/ferroptosis", "Fibrotic", "Composite")
    )
  )
write_csv(module_trend_df, file.path(table_dir, "npc_slingshot_module_trends_along_pseudotime.csv"))

p_module_trends <- ggplot(module_trend_df, aes(pseudotime, mean_score, color = module)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("#4E79A7", "#F28E2B", "#E15759", "#6B6B6B")) +
  labs(
    title = "F. Module trends along Slingshot pseudotime",
    subtitle = "Binned mean scaled module scores; descriptive trend only",
    x = "Slingshot pseudotime",
    y = "Mean module score",
    color = NULL
  ) +
  theme_pt(9) +
  theme(legend.position = "none")

main_fig <- (p_lineage | p_pt) / (p_box_global | p_module_trends) +
  plot_annotation(
    title = "Slingshot analysis supports a transcriptional continuum across NPC states",
    subtitle = "Trajectory inference is exploratory and does not establish causal or real-time progression",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )

ggsave(file.path(figure_dir, "51_npc_slingshot_embedding_lineages.png"), p_lineage, width = 7.2, height = 5.8, dpi = 300)
ggsave(file.path(figure_dir, "52_npc_slingshot_pseudotime_embedding.png"), p_pt, width = 7.2, height = 5.8, dpi = 300)
ggsave(file.path(figure_dir, "53_npc_slingshot_pseudotime_by_group.png"), p_box_global, width = 5.2, height = 4.6, dpi = 300)
ggsave(file.path(figure_dir, "54_npc_slingshot_pseudotime_by_subtype_group.png"), p_box_subtype, width = 11, height = 5.6, dpi = 300)
ggsave(file.path(figure_dir, "55_npc_slingshot_module_trends.png"), p_module_trends, width = 12, height = 4.0, dpi = 300)
ggsave(file.path(figure_dir, "56_npc_slingshot_main_Figure6_candidate.png"), main_fig, width = 13.5, height = 9.2, dpi = 300)

saveRDS(sce, file.path(processed_dir, "GSE165722_clean_npc_subset_slingshot_sce.rds"))

message("Slingshot pseudotime analysis finished.")
