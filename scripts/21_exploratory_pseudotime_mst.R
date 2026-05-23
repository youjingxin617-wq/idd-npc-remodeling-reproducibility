suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(igraph)
  library(ggrepel)
  library(pheatmap)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures_pseudotime")
table_dir <- file.path(project_root, "results/tables_pseudotime")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260521)

pal_group <- c(mild = "#4E79A7", severe = "#E15759")
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

theme_pt <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

get_data_layer <- function(obj, features = NULL) {
  mat <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
  )
  if (!is.null(features)) {
    features <- intersect(features, rownames(mat))
    mat <- mat[features, , drop = FALSE]
  }
  mat
}

project_to_segment <- function(x, a, b) {
  ab <- b - a
  denom <- sum(ab * ab)
  if (denom == 0) {
    t <- rep(0, nrow(x))
  } else {
    t <- ((x[, 1] - a[1]) * ab[1] + (x[, 2] - a[2]) * ab[2]) / denom
    t <- pmax(0, pmin(1, t))
  }
  proj <- cbind(a[1] + t * ab[1], a[2] + t * ab[2])
  dist <- sqrt(rowSums((x - proj)^2))
  list(t = t, dist = dist)
}

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"

reduction_use <- if ("umap" %in% Reductions(clean_npc)) "umap" else "tsne"
emb <- as.data.frame(Embeddings(clean_npc, reduction = reduction_use)[, 1:2, drop = FALSE])
colnames(emb) <- c("dim1", "dim2")

meta <- clean_npc@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  bind_cols(emb) %>%
  mutate(
    npc_subtype_final = factor(npc_subtype_final, levels = subtype_levels),
    article_group = factor(article_group, levels = c("mild", "severe"))
  ) %>%
  filter(!is.na(npc_subtype_final), !is.na(article_group))

centroids <- meta %>%
  group_by(npc_subtype_final) %>%
  summarise(
    dim1 = median(dim1, na.rm = TRUE),
    dim2 = median(dim2, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(npc_subtype_final)

centroid_mat <- as.matrix(centroids[, c("dim1", "dim2")])
rownames(centroid_mat) <- as.character(centroids$npc_subtype_final)
dist_mat <- as.matrix(dist(centroid_mat))
g <- graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
mst_g <- mst(g, weights = E(g)$weight)

root_subtype <- "Hom/Reg-like NPCs"
if (!root_subtype %in% V(mst_g)$name) {
  root_subtype <- as.character(centroids$npc_subtype_final[which.max(centroids$n_cells)])
}
root_dist <- distances(mst_g, v = root_subtype, weights = E(mst_g)$weight)[1, ]

edge_tbl <- as_data_frame(mst_g, what = "edges") %>%
  mutate(
    from_dist = as.numeric(root_dist[from]),
    to_dist = as.numeric(root_dist[to]),
    parent = if_else(from_dist <= to_dist, from, to),
    child = if_else(from_dist <= to_dist, to, from),
    parent_dist = pmin(from_dist, to_dist),
    child_dist = pmax(from_dist, to_dist)
  ) %>%
  rowwise() %>%
  mutate(
    parent_dim1 = centroid_mat[parent, "dim1"],
    parent_dim2 = centroid_mat[parent, "dim2"],
    child_dim1 = centroid_mat[child, "dim1"],
    child_dim2 = centroid_mat[child, "dim2"],
    edge_length = sqrt((child_dim1 - parent_dim1)^2 + (child_dim2 - parent_dim2)^2)
  ) %>%
  ungroup()

x <- as.matrix(meta[, c("dim1", "dim2")])
best_dist <- rep(Inf, nrow(meta))
best_pt <- rep(NA_real_, nrow(meta))
best_edge <- rep(NA_character_, nrow(meta))

for (i in seq_len(nrow(edge_tbl))) {
  a <- c(edge_tbl$parent_dim1[i], edge_tbl$parent_dim2[i])
  b <- c(edge_tbl$child_dim1[i], edge_tbl$child_dim2[i])
  proj <- project_to_segment(x, a, b)
  pt <- edge_tbl$parent_dist[i] + proj$t * edge_tbl$edge_length[i]
  replace <- proj$dist < best_dist
  best_dist[replace] <- proj$dist[replace]
  best_pt[replace] <- pt[replace]
  best_edge[replace] <- paste(edge_tbl$parent[i], edge_tbl$child[i], sep = " -> ")
}

meta$pseudotime_raw <- best_pt
meta$pseudotime <- (best_pt - min(best_pt, na.rm = TRUE)) / diff(range(best_pt, na.rm = TRUE))
meta$assigned_mst_edge <- best_edge
meta$distance_to_mst <- best_dist

write_csv(meta, file.path(table_dir, "npc_mst_exploratory_pseudotime_per_cell.csv"))
write_csv(centroids, file.path(table_dir, "npc_mst_subtype_centroids.csv"))
write_csv(edge_tbl, file.path(table_dir, "npc_mst_edges_rooted_from_Hom_Reg_like_NPCs.csv"))

pt_summary_subtype <- meta %>%
  group_by(npc_subtype_final, article_group) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    sd_pseudotime = sd(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(pt_summary_subtype, file.path(table_dir, "npc_mst_pseudotime_summary_by_subtype_group.csv"))

pt_stats_subtype <- meta %>%
  group_by(npc_subtype_final) %>%
  summarise(
    n_mild = sum(article_group == "mild"),
    n_severe = sum(article_group == "severe"),
    mean_mild = mean(pseudotime[article_group == "mild"], na.rm = TRUE),
    mean_severe = mean(pseudotime[article_group == "severe"], na.rm = TRUE),
    delta_severe_minus_mild = mean_severe - mean_mild,
    p_value = {
      mild_values <- pseudotime[article_group == "mild"]
      severe_values <- pseudotime[article_group == "severe"]
      if (length(mild_values) >= 3 && length(severe_values) >= 3) {
        wilcox.test(mild_values, severe_values, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
write_csv(pt_stats_subtype, file.path(table_dir, "npc_mst_pseudotime_mild_vs_severe_stats_by_subtype.csv"))

pt_stats_global <- meta %>%
  summarise(
    n_mild = sum(article_group == "mild"),
    n_severe = sum(article_group == "severe"),
    mean_mild = mean(pseudotime[article_group == "mild"], na.rm = TRUE),
    mean_severe = mean(pseudotime[article_group == "severe"], na.rm = TRUE),
    delta_severe_minus_mild = mean_severe - mean_mild,
    p_value = wilcox.test(
      pseudotime[article_group == "mild"],
      pseudotime[article_group == "severe"],
      exact = FALSE
    )$p.value
  )
write_csv(pt_stats_global, file.path(table_dir, "npc_mst_pseudotime_global_mild_vs_severe_stats.csv"))

line_df <- edge_tbl %>%
  transmute(
    x = parent_dim1,
    y = parent_dim2,
    xend = child_dim1,
    yend = child_dim2,
    parent,
    child
  )

p_subtype <- ggplot(meta, aes(dim1, dim2, color = npc_subtype_final)) +
  geom_point(size = 0.16, alpha = 0.48) +
  geom_segment(data = line_df, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, color = "black", linewidth = 0.65, arrow = arrow(length = unit(0.10, "inches"))) +
  geom_point(data = centroids, aes(dim1, dim2), inherit.aes = FALSE, size = 2.4, color = "black", fill = "white", shape = 21, stroke = 0.7) +
  geom_text_repel(data = centroids, aes(dim1, dim2, label = npc_subtype_final), inherit.aes = FALSE, size = 2.8, max.overlaps = 20) +
  scale_color_manual(values = subtype_pal, drop = FALSE) +
  labs(
    title = "A. MST-based exploratory trajectory across NPC states",
    subtitle = paste0("Root set to ", root_subtype, "; arrows show centroid-level MST direction"),
    x = paste0(toupper(reduction_use), " 1"),
    y = paste0(toupper(reduction_use), " 2"),
    color = "NPC subtype"
  ) +
  theme_pt(9)

p_pt <- ggplot(meta, aes(dim1, dim2, color = pseudotime)) +
  geom_point(size = 0.16, alpha = 0.62) +
  geom_segment(data = line_df, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, color = "black", linewidth = 0.55, arrow = arrow(length = unit(0.10, "inches"))) +
  scale_color_gradientn(colors = c("#2D5F9A", "#F4E8B5", "#B53636")) +
  labs(
    title = "B. Continuous exploratory pseudotime",
    subtitle = "Each NPC is projected to the nearest MST segment",
    x = paste0(toupper(reduction_use), " 1"),
    y = paste0(toupper(reduction_use), " 2"),
    color = "Pseudotime"
  ) +
  theme_pt(9)

p_group <- ggplot(meta, aes(dim1, dim2, color = article_group)) +
  geom_point(size = 0.16, alpha = 0.52) +
  geom_segment(data = line_df, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, color = "black", linewidth = 0.55, arrow = arrow(length = unit(0.10, "inches"))) +
  scale_color_manual(values = pal_group) +
  labs(
    title = "C. Degeneration group along exploratory trajectory",
    subtitle = "Distribution is descriptive and should not be interpreted as real temporal progression",
    x = paste0(toupper(reduction_use), " 1"),
    y = paste0(toupper(reduction_use), " 2"),
    color = "Group"
  ) +
  theme_pt(9)

combined_embedding <- (p_subtype | p_pt | p_group) +
  plot_annotation(
    title = "Exploratory MST-based pseudotime analysis of NPC states",
    subtitle = "Trajectory inference is hypothesis-generating and does not establish true temporal or causal order",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )
ggsave(file.path(figure_dir, "44_npc_mst_pseudotime_embedding_combined.png"), combined_embedding, width = 14.5, height = 5.0, dpi = 300)

p_box_subtype <- ggplot(meta, aes(npc_subtype_final, pseudotime, fill = article_group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.82, linewidth = 0.3, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = article_group), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 0.18, alpha = 0.18, show.legend = FALSE) +
  scale_fill_manual(values = pal_group) +
  scale_color_manual(values = pal_group) +
  labs(
    title = "Exploratory pseudotime distribution by NPC subtype",
    subtitle = "Root: Hom/Reg-like NPCs; higher values indicate later positions along the inferred MST",
    x = NULL,
    y = "MST-based pseudotime",
    fill = "Group"
  ) +
  theme_pt(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(figure_dir, "45_npc_mst_pseudotime_by_subtype_group.png"), p_box_subtype, width = 11, height = 5.6, dpi = 300)

p_box_global <- ggplot(meta, aes(article_group, pseudotime, fill = article_group)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.82, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.18, alpha = 0.16) +
  scale_fill_manual(values = pal_group) +
  labs(
    title = "Global exploratory pseudotime by degeneration group",
    subtitle = "Cell-level comparison; interpret together with subtype-composition changes",
    x = NULL,
    y = "MST-based pseudotime",
    fill = "Group"
  ) +
  theme_pt(10) +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "46_npc_mst_pseudotime_global_mild_vs_severe.png"), p_box_global, width = 5.2, height = 4.6, dpi = 300)

core_gene_sets <- list(
  Matrix_homeostasis = c("ACAN", "COL2A1", "COL6A2", "BGN"),
  Iron_ferroptosis_response = c("GPX4", "SLC39A14", "FTH1", "FTL", "SLC7A11"),
  Fibrotic_remodeling = c("COL1A1", "COL3A1", "POSTN", "ADAMTS5", "MMP13")
)
trend_genes <- c("ACAN", "COL2A1", "COL1A1", "COL3A1", "POSTN", "MMP13", "GPX4", "SLC39A14", "FTH1", "FTL")
genes_use <- intersect(unique(c(trend_genes, unlist(core_gene_sets))), rownames(clean_npc))

expr <- as.matrix(get_data_layer(clean_npc, genes_use)[, meta$cell_barcode, drop = FALSE])
scaled_expr <- t(scale(t(expr)))
scaled_expr[is.na(scaled_expr)] <- 0

trend_df <- as.data.frame(t(scaled_expr[intersect(trend_genes, rownames(scaled_expr)), , drop = FALSE])) %>%
  rownames_to_column("cell_barcode") %>%
  left_join(meta %>% select(cell_barcode, pseudotime, article_group, npc_subtype_final), by = "cell_barcode") %>%
  mutate(pseudotime_bin = cut(pseudotime, breaks = seq(0, 1, length.out = 26), include.lowest = TRUE, labels = FALSE)) %>%
  pivot_longer(all_of(intersect(trend_genes, rownames(scaled_expr))), names_to = "gene", values_to = "z_expr") %>%
  group_by(gene, pseudotime_bin) %>%
  summarise(
    pseudotime = mean(pseudotime, na.rm = TRUE),
    mean_z = mean(z_expr, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(trend_df, file.path(table_dir, "npc_mst_core_gene_binned_expression_along_pseudotime.csv"))

p_gene_trends <- ggplot(trend_df, aes(pseudotime, mean_z, color = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_color_manual(values = rep(c("#4E79A7", "#E15759", "#59A14F", "#F28E2B", "#B07AA1", "#9C755F"), length.out = length(unique(trend_df$gene)))) +
  labs(
    title = "Core gene expression trends along exploratory pseudotime",
    subtitle = "Binned mean z-scored expression; descriptive trend only",
    x = "MST-based pseudotime",
    y = "Mean scaled expression",
    color = "Gene"
  ) +
  theme_pt(9) +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "47_npc_mst_core_gene_trends_along_pseudotime.png"), p_gene_trends, width = 12, height = 6.2, dpi = 300)

module_mat <- map_dfc(names(core_gene_sets), function(module_name) {
  genes <- intersect(core_gene_sets[[module_name]], rownames(scaled_expr))
  tibble(!!module_name := colMeans(scaled_expr[genes, , drop = FALSE]))
})
module_df <- tibble(cell_barcode = colnames(scaled_expr)) %>%
  bind_cols(module_mat) %>%
  left_join(meta %>% select(cell_barcode, pseudotime, article_group, npc_subtype_final), by = "cell_barcode") %>%
  mutate(
    Composite_remodeling_iron_over_matrix = Iron_ferroptosis_response + Fibrotic_remodeling - Matrix_homeostasis,
    pseudotime_bin = cut(pseudotime, breaks = seq(0, 1, length.out = 26), include.lowest = TRUE, labels = FALSE)
  )
write_csv(module_df, file.path(table_dir, "npc_mst_module_scores_per_cell.csv"))

module_trend_df <- module_df %>%
  pivot_longer(c(Matrix_homeostasis, Iron_ferroptosis_response, Fibrotic_remodeling, Composite_remodeling_iron_over_matrix), names_to = "module", values_to = "score") %>%
  group_by(module, pseudotime_bin) %>%
  summarise(
    pseudotime = mean(pseudotime, na.rm = TRUE),
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
write_csv(module_trend_df, file.path(table_dir, "npc_mst_module_trends_along_pseudotime.csv"))

p_module_trends <- ggplot(module_trend_df, aes(pseudotime, mean_score, color = module)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("#4E79A7", "#F28E2B", "#E15759", "#6B6B6B")) +
  labs(
    title = "Core module trends along exploratory pseudotime",
    subtitle = "Binned mean scaled module scores; smoothed for visualization",
    x = "MST-based pseudotime",
    y = "Mean module score",
    color = NULL
  ) +
  theme_pt(9) +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "48_npc_mst_module_trends_along_pseudotime.png"), p_module_trends, width = 12, height = 4.0, dpi = 300)

heat_df <- trend_df %>%
  select(gene, pseudotime_bin, mean_z) %>%
  pivot_wider(names_from = pseudotime_bin, values_from = mean_z) %>%
  arrange(match(gene, trend_genes))
heat_mat <- heat_df %>% column_to_rownames("gene") %>% as.matrix()
png(file.path(figure_dir, "49_npc_mst_core_gene_pseudotime_heatmap.png"), width = 1800, height = 900, res = 200)
pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2D5F9A", "white", "#B53636"))(100),
  fontsize = 8,
  main = "Core gene trends along MST-based pseudotime"
)
dev.off()

main_fig <- (p_subtype | p_pt) / (p_box_global | p_module_trends) +
  plot_annotation(
    title = "Exploratory pseudotime suggests a transcriptional continuum across NPC states",
    subtitle = "MST-based analysis rooted at Hom/Reg-like NPCs; results are hypothesis-generating",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )
ggsave(file.path(figure_dir, "50_npc_mst_pseudotime_main_Figure6_candidate.png"), main_fig, width = 13.5, height = 9.2, dpi = 300)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("Exploratory MST-based pseudotime analysis finished.")
