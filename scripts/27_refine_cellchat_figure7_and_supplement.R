suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

project_root <- normalizePath(".", winslash = "/")
table_dir <- file.path(project_root, "results/tables_cellchat")
figure_dir <- file.path(project_root, "results/figures_cellchat")
vault_figure_dir <- "E:/github/youwenblood/生信复现sci文献优化/06_细胞通讯分析/figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vault_figure_dir, recursive = TRUE, showWarnings = FALSE)

subtype_levels <- c(
  "ECM/Adh-NPCs", "Hom-NPCs", "Eff-NPCs", "Ht-NPCs",
  "Fibro-NPCs", "Fibro-reg NPCs", "Hom/Reg-like NPCs"
)
pal_group <- c(mild = "#4E79A7", severe = "#E15759")

theme_comm <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

net_long <- read_csv(file.path(table_dir, "npc_cellchat_network_count_weight_matrices_long.csv"), show_col_types = FALSE) %>%
  mutate(
    source = factor(source, levels = subtype_levels),
    target = factor(target, levels = subtype_levels)
  )

centrality_tbl <- read_csv(file.path(table_dir, "npc_cellchat_outgoing_incoming_summary_by_subtype.csv"), show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = c("mild", "severe")),
    subtype = factor(subtype, levels = subtype_levels)
  )

top_pathways <- read_csv(file.path(table_dir, "npc_cellchat_top_severe_enhanced_pathways.csv"), show_col_types = FALSE) %>%
  slice_head(n = 12) %>%
  mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))

top_pairs <- read_csv(file.path(table_dir, "npc_cellchat_top_severe_enhanced_lr_pairs.csv"), show_col_types = FALSE) %>%
  slice_head(n = 24) %>%
  mutate(
    source = factor(source, levels = rev(subtype_levels)),
    target = factor(target, levels = subtype_levels),
    pair_label = paste0(interaction_name_2, "\n", pathway_name),
    pair_label = factor(pair_label, levels = rev(unique(pair_label)))
  )

p_delta_heat <- net_long %>%
  filter(metric == "weight_delta", group == "severe_minus_mild") %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2D5F9A", mid = "white", high = "#B53636", midpoint = 0) +
  labs(
    title = "A. Severe - mild inferred communication strength",
    subtitle = "Positive values indicate stronger inferred signaling in severe NPCs",
    x = "Target subtype",
    y = "Source subtype",
    fill = "Delta"
  ) +
  theme_comm(8)

centrality_plot_df <- centrality_tbl %>%
  select(group, subtype, outgoing_weight, incoming_weight) %>%
  pivot_longer(c(outgoing_weight, incoming_weight), names_to = "direction", values_to = "weight") %>%
  mutate(direction = recode(direction, outgoing_weight = "Outgoing", incoming_weight = "Incoming"))

p_centrality <- ggplot(centrality_plot_df, aes(subtype, weight, fill = group)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  facet_wrap(~ direction, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_group) +
  labs(
    title = "B. NPC subtype incoming and outgoing strength",
    x = NULL,
    y = "Aggregated CellChat weight",
    fill = "Group"
  ) +
  theme_comm(8)

p_pathways <- ggplot(top_pathways, aes(total_prob_delta_severe_minus_mild, pathway_name)) +
  geom_col(fill = "#B53636", alpha = 0.84) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#666666", linewidth = 0.3) +
  labs(
    title = "C. Severe-enhanced inferred pathways",
    x = "Delta total communication probability",
    y = NULL
  ) +
  theme_comm(8) +
  theme(axis.text.x = element_text(angle = 0))

main_fig <- p_delta_heat | (p_centrality / p_pathways) +
  plot_layout(widths = c(1.1, 1.25)) +
  plot_annotation(
    title = "Exploratory CellChat analysis suggests ECM-centered NPC communication remodeling",
    subtitle = "Database-dependent ligand-receptor inference from scRNA-seq; not experimental proof of ligand-receptor activity",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )

p_count_heat <- net_long %>%
  filter(metric == "count", group %in% c("mild", "severe")) %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
  labs(title = "A. Number of inferred communications", x = "Target", y = "Source", fill = "Count") +
  theme_comm(8)

p_weight_heat <- net_long %>%
  filter(metric == "weight", group %in% c("mild", "severe")) %>%
  ggplot(aes(target, source, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_gradient(low = "#FFF5EB", high = "#A63603") +
  labs(title = "B. Inferred communication strength", x = "Target", y = "Source", fill = "Weight") +
  theme_comm(8)

p_pairs <- ggplot(top_pairs, aes(x = target, y = pair_label)) +
  geom_point(aes(size = prob_severe, color = prob_delta_severe_minus_mild), alpha = 0.86) +
  facet_wrap(~ source, scales = "free_y") +
  scale_color_gradient2(low = "#2D5F9A", mid = "#F4E8B5", high = "#B53636", midpoint = 0) +
  labs(
    title = "C. Selected severe-enhanced ligand-receptor interactions",
    subtitle = "Ranked by severe-minus-mild communication probability",
    x = "Target subtype",
    y = "Ligand-receptor pair / pathway",
    size = "Severe prob.",
    color = "Delta prob."
  ) +
  theme_comm(7)

supp_fig <- (p_count_heat / p_weight_heat / p_pairs) +
  plot_layout(heights = c(1, 1, 1.65)) +
  plot_annotation(
    title = "Supplementary CellChat outputs supporting Figure 7",
    subtitle = "Full count and weight heatmaps plus selected severe-enhanced ligand-receptor pairs",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  )

out_main <- file.path(figure_dir, "63_npc_cellchat_main_Figure7_candidate.png")
out_supp <- file.path(figure_dir, "64_npc_cellchat_supplementary_FigureS1_candidate.png")
ggsave(out_main, main_fig, width = 13.5, height = 8.0, dpi = 300)
ggsave(out_supp, supp_fig, width = 14.0, height = 16.0, dpi = 300)

file.copy(out_main, file.path(vault_figure_dir, basename(out_main)), overwrite = TRUE)
file.copy(out_supp, file.path(vault_figure_dir, basename(out_supp)), overwrite = TRUE)

message("Wrote refined Figure 7: ", out_main)
message("Wrote Supplementary Figure S1: ", out_supp)
