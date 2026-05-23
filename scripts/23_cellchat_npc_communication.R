suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
  library(igraph)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures_cellchat")
table_dir <- file.path(project_root, "results/tables_cellchat")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260521)
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

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

theme_comm <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

mat_to_long <- function(mat, group_name, metric) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    set_names(c("source", "target", "value")) %>%
    mutate(group = group_name, metric = metric)
}

run_cellchat_one_group <- function(seurat_obj, group_name, db_use, min_cells = 10) {
  message("Running CellChat for ", group_name, " ...")
  obj_group <- subset(seurat_obj, subset = article_group == group_name)
  obj_group$cellchat_label <- factor(obj_group$npc_subtype_final, levels = subtype_levels)
  data.input <- safe_get_assay(obj_group, layer = "data")
  meta <- obj_group[[]] %>%
    transmute(labels = cellchat_label) %>%
    as.data.frame()
  rownames(meta) <- colnames(obj_group)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "RNA")
  cellchat@DB <- db_use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = TRUE,
    nboot = 20,
    seed.use = 20260521
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"
clean_npc$npc_subtype_final <- factor(clean_npc$npc_subtype_final, levels = subtype_levels)
clean_npc$article_group <- factor(clean_npc$article_group, levels = c("mild", "severe"))
clean_npc <- subset(clean_npc, subset = !is.na(npc_subtype_final) & !is.na(article_group))

cell_count_tbl <- clean_npc[[]] %>%
  count(article_group, npc_subtype_final, name = "n_cells") %>%
  complete(article_group, npc_subtype_final = factor(subtype_levels, levels = subtype_levels), fill = list(n_cells = 0))
write_csv(cell_count_tbl, file.path(table_dir, "npc_cellchat_cell_counts_by_group_subtype.csv"))

data(CellChatDB.human)
# Full human database is used because NPC remodeling may involve secreted signaling, ECM-receptor,
# and cell-cell contact interactions. This stays exploratory and database-dependent.
CellChatDB.use <- CellChatDB.human

cellchat_mild <- run_cellchat_one_group(clean_npc, "mild", CellChatDB.use)
cellchat_severe <- run_cellchat_one_group(clean_npc, "severe", CellChatDB.use)

saveRDS(cellchat_mild, file.path(processed_dir, "GSE165722_clean_npc_CellChat_mild.rds"))
saveRDS(cellchat_severe, file.path(processed_dir, "GSE165722_clean_npc_CellChat_severe.rds"))

comm_mild <- subsetCommunication(cellchat_mild) %>% mutate(group = "mild")
comm_severe <- subsetCommunication(cellchat_severe) %>% mutate(group = "severe")
comm_all <- bind_rows(comm_mild, comm_severe)
write_csv(comm_all, file.path(table_dir, "npc_cellchat_all_significant_interactions_mild_severe.csv"))

path_mild <- subsetCommunication(cellchat_mild, slot.name = "netP") %>% mutate(group = "mild")
path_severe <- subsetCommunication(cellchat_severe, slot.name = "netP") %>% mutate(group = "severe")
path_all <- bind_rows(path_mild, path_severe)
write_csv(path_all, file.path(table_dir, "npc_cellchat_significant_pathway_interactions_mild_severe.csv"))

count_mild <- cellchat_mild@net$count
count_severe <- cellchat_severe@net$count
weight_mild <- cellchat_mild@net$weight
weight_severe <- cellchat_severe@net$weight

# Align matrices to fixed subtype order.
count_mild <- count_mild[subtype_levels, subtype_levels, drop = FALSE]
count_severe <- count_severe[subtype_levels, subtype_levels, drop = FALSE]
weight_mild <- weight_mild[subtype_levels, subtype_levels, drop = FALSE]
weight_severe <- weight_severe[subtype_levels, subtype_levels, drop = FALSE]

net_long <- bind_rows(
  mat_to_long(count_mild, "mild", "count"),
  mat_to_long(count_severe, "severe", "count"),
  mat_to_long(count_severe - count_mild, "severe_minus_mild", "count_delta"),
  mat_to_long(weight_mild, "mild", "weight"),
  mat_to_long(weight_severe, "severe", "weight"),
  mat_to_long(weight_severe - weight_mild, "severe_minus_mild", "weight_delta")
) %>%
  mutate(
    source = factor(source, levels = subtype_levels),
    target = factor(target, levels = subtype_levels)
  )
write_csv(net_long, file.path(table_dir, "npc_cellchat_network_count_weight_matrices_long.csv"))

summary_total <- tibble(
  group = c("mild", "severe"),
  total_interaction_count = c(sum(count_mild), sum(count_severe)),
  total_interaction_weight = c(sum(weight_mild), sum(weight_severe)),
  mean_interaction_weight = c(mean(weight_mild[weight_mild > 0]), mean(weight_severe[weight_severe > 0]))
)
write_csv(summary_total, file.path(table_dir, "npc_cellchat_global_network_summary.csv"))

centrality_tbl <- bind_rows(
  tibble(group = "mild", subtype = subtype_levels, outgoing_weight = rowSums(weight_mild), incoming_weight = colSums(weight_mild), outgoing_count = rowSums(count_mild), incoming_count = colSums(count_mild)),
  tibble(group = "severe", subtype = subtype_levels, outgoing_weight = rowSums(weight_severe), incoming_weight = colSums(weight_severe), outgoing_count = rowSums(count_severe), incoming_count = colSums(count_severe))
) %>%
  mutate(subtype = factor(subtype, levels = subtype_levels))
write_csv(centrality_tbl, file.path(table_dir, "npc_cellchat_outgoing_incoming_summary_by_subtype.csv"))

centrality_delta <- centrality_tbl %>%
  pivot_wider(names_from = group, values_from = c(outgoing_weight, incoming_weight, outgoing_count, incoming_count)) %>%
  mutate(
    outgoing_weight_delta = outgoing_weight_severe - outgoing_weight_mild,
    incoming_weight_delta = incoming_weight_severe - incoming_weight_mild,
    outgoing_count_delta = outgoing_count_severe - outgoing_count_mild,
    incoming_count_delta = incoming_count_severe - incoming_count_mild
  )
write_csv(centrality_delta, file.path(table_dir, "npc_cellchat_outgoing_incoming_delta_by_subtype.csv"))

pair_compare <- comm_all %>%
  mutate(pair_key = paste(source, target, interaction_name_2, pathway_name, sep = "|")) %>%
  group_by(group, source, target, interaction_name_2, pathway_name, pair_key) %>%
  summarise(prob = max(prob, na.rm = TRUE), pval = min(pval, na.rm = TRUE), .groups = "drop") %>%
  select(group, source, target, interaction_name_2, pathway_name, pair_key, prob, pval) %>%
  pivot_wider(names_from = group, values_from = c(prob, pval), values_fill = list(prob = 0, pval = NA_real_)) %>%
  mutate(
    prob_delta_severe_minus_mild = prob_severe - prob_mild,
    severe_specific_or_enhanced = prob_severe > 0 & prob_delta_severe_minus_mild > 0
  ) %>%
  arrange(desc(prob_delta_severe_minus_mild), desc(prob_severe))
write_csv(pair_compare, file.path(table_dir, "npc_cellchat_ligand_receptor_pair_comparison_mild_vs_severe.csv"))

top_pairs <- pair_compare %>%
  filter(severe_specific_or_enhanced) %>%
  slice_head(n = 35) %>%
  mutate(
    source = factor(source, levels = rev(subtype_levels)),
    target = factor(target, levels = subtype_levels),
    pair_label = paste0(interaction_name_2, "\n", pathway_name),
    pair_label = factor(pair_label, levels = rev(unique(pair_label)))
  )
write_csv(top_pairs, file.path(table_dir, "npc_cellchat_top_severe_enhanced_lr_pairs.csv"))

pathway_compare <- path_all %>%
  group_by(group, pathway_name) %>%
  summarise(n_pairs = n(), total_prob = sum(prob, na.rm = TRUE), mean_prob = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(n_pairs, total_prob, mean_prob), values_fill = 0) %>%
  mutate(
    total_prob_delta_severe_minus_mild = total_prob_severe - total_prob_mild,
    n_pairs_delta_severe_minus_mild = n_pairs_severe - n_pairs_mild
  ) %>%
  arrange(desc(total_prob_delta_severe_minus_mild))
write_csv(pathway_compare, file.path(table_dir, "npc_cellchat_pathway_comparison_mild_vs_severe.csv"))

top_pathways <- pathway_compare %>%
  slice_head(n = 25) %>%
  mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))
write_csv(top_pathways, file.path(table_dir, "npc_cellchat_top_severe_enhanced_pathways.csv"))

p_count_heat <- net_long %>%
  filter(metric %in% c("count"), group %in% c("mild", "severe")) %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
  labs(title = "A. Number of inferred NPC subtype communications", x = "Target", y = "Source", fill = "Count") +
  theme_comm(8)

p_weight_heat <- net_long %>%
  filter(metric %in% c("weight"), group %in% c("mild", "severe")) %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradient(low = "#FFF5EB", high = "#A63603") +
  labs(title = "B. Inferred communication strength", x = "Target", y = "Source", fill = "Weight") +
  theme_comm(8)

p_delta_heat <- net_long %>%
  filter(metric == "weight_delta", group == "severe_minus_mild") %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2D5F9A", mid = "white", high = "#B53636", midpoint = 0) +
  labs(title = "C. Severe - mild communication strength", subtitle = "Positive values indicate stronger inferred signaling in severe NPCs", x = "Target", y = "Source", fill = "Delta") +
  theme_comm(8)

centrality_plot_df <- centrality_tbl %>%
  select(group, subtype, outgoing_weight, incoming_weight) %>%
  pivot_longer(c(outgoing_weight, incoming_weight), names_to = "direction", values_to = "weight") %>%
  mutate(direction = recode(direction, outgoing_weight = "Outgoing", incoming_weight = "Incoming"))

p_centrality <- ggplot(centrality_plot_df, aes(subtype, weight, fill = group)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  facet_wrap(~ direction, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_group) +
  labs(title = "D. NPC subtype outgoing and incoming communication strength", x = NULL, y = "Aggregated CellChat weight", fill = "Group") +
  theme_comm(8)

p_pairs <- ggplot(top_pairs, aes(x = target, y = pair_label)) +
  geom_point(aes(size = prob_severe, color = prob_delta_severe_minus_mild), alpha = 0.85) +
  facet_wrap(~ source, scales = "free_y") +
  scale_color_gradient2(low = "#2D5F9A", mid = "#F4E8B5", high = "#B53636", midpoint = 0) +
  labs(
    title = "E. Top severe-enhanced ligand-receptor interactions",
    subtitle = "Selected by severe-minus-mild communication probability; exploratory CellChat inference",
    x = "Target subtype",
    y = "Ligand-receptor pair / pathway",
    size = "Severe prob.",
    color = "Delta prob."
  ) +
  theme_comm(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pathways <- ggplot(top_pathways, aes(total_prob_delta_severe_minus_mild, pathway_name)) +
  geom_col(fill = "#B53636", alpha = 0.82) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#666666", linewidth = 0.3) +
  labs(
    title = "F. Pathways with stronger inferred communication in severe NPCs",
    x = "Delta total communication probability",
    y = NULL
  ) +
  theme_comm(8) +
  theme(axis.text.x = element_text(angle = 0))

combined_main <- (p_count_heat / p_delta_heat) | (p_centrality / p_pathways) +
  plot_annotation(
    title = "Exploratory CellChat analysis of communication among NPC subtypes",
    subtitle = "Database-dependent inference from scRNA-seq expression; not experimental proof of ligand-receptor activity",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )

ggsave(file.path(figure_dir, "57_npc_cellchat_count_heatmaps.png"), p_count_heat, width = 10.5, height = 4.5, dpi = 300)
ggsave(file.path(figure_dir, "58_npc_cellchat_weight_heatmaps.png"), p_weight_heat, width = 10.5, height = 4.5, dpi = 300)
ggsave(file.path(figure_dir, "59_npc_cellchat_weight_delta_heatmap.png"), p_delta_heat, width = 6.0, height = 5.2, dpi = 300)
ggsave(file.path(figure_dir, "60_npc_cellchat_outgoing_incoming_strength.png"), p_centrality, width = 10.5, height = 4.8, dpi = 300)
ggsave(file.path(figure_dir, "61_npc_cellchat_top_severe_enhanced_lr_pairs.png"), p_pairs, width = 14.5, height = 9.5, dpi = 300)
ggsave(file.path(figure_dir, "62_npc_cellchat_top_severe_enhanced_pathways.png"), p_pathways, width = 7.0, height = 6.2, dpi = 300)
ggsave(file.path(figure_dir, "63_npc_cellchat_main_Figure7_candidate.png"), combined_main, width = 14.0, height = 10.0, dpi = 300)

message("CellChat NPC communication analysis finished.")
