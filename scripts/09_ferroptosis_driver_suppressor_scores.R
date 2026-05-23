suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures")
table_dir <- file.path(project_root, "results/tables")
gene_set_dir <- file.path(project_root, "data/gene_sets/FerrDb")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260519)

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"

driver_file <- file.path(gene_set_dir, "ferroptosis_driver.csv")
suppressor_file <- file.path(gene_set_dir, "ferroptosis_suppressor.csv")

if (!file.exists(driver_file) || !file.exists(suppressor_file)) {
  stop(
    "FerrDb gene-set files were not found. Expected:\n",
    driver_file, "\n",
    suppressor_file, "\n",
    "Please download FerrDb ferroptosis driver/suppressor CSV files first."
  )
}

read_ferrdb_genes <- function(path) {
  ferrdb <- readr::read_csv(path, show_col_types = FALSE)
  ferrdb %>%
    filter(!is.na(symbol), symbol != "") %>%
    filter(is.na(testin) | str_detect(testin, regex("human", ignore_case = TRUE))) %>%
    pull(symbol) %>%
    unique()
}

driver_genes_all <- read_ferrdb_genes(driver_file)
suppressor_genes_all <- read_ferrdb_genes(suppressor_file)

driver_genes <- intersect(driver_genes_all, rownames(clean_npc))
suppressor_genes <- intersect(suppressor_genes_all, rownames(clean_npc))

if (length(driver_genes) < 5) {
  stop("Too few FerrDb driver genes were detected in the expression matrix: ", length(driver_genes))
}
if (length(suppressor_genes) < 5) {
  stop("Too few FerrDb suppressor genes were detected in the expression matrix: ", length(suppressor_genes))
}

message("FerrDb driver genes downloaded/listed: ", length(driver_genes_all))
message("FerrDb driver genes found in this NPC matrix: ", length(driver_genes))
message("FerrDb suppressor genes downloaded/listed: ", length(suppressor_genes_all))
message("FerrDb suppressor genes found in this NPC matrix: ", length(suppressor_genes))

clean_npc <- AddModuleScore(
  clean_npc,
  features = list(driver_genes),
  assay = "RNA",
  name = "FDS"
)
clean_npc$FDS <- clean_npc$FDS1

clean_npc <- AddModuleScore(
  clean_npc,
  features = list(suppressor_genes),
  assay = "RNA",
  name = "FSS"
)
clean_npc$FSS <- clean_npc$FSS1

score_df <- clean_npc@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  select(
    cell_barcode,
    sample_id,
    grade,
    article_group,
    npc_clean_cluster,
    npc_subtype_final,
    FDS,
    FSS
  ) %>%
  mutate(
    article_group = factor(article_group, levels = c("mild", "severe")),
    npc_subtype_final = factor(
      npc_subtype_final,
      levels = c(
        "ECM/Adh-NPCs",
        "Hom-NPCs",
        "Eff-NPCs",
        "Ht-NPCs",
        "Fibro-NPCs",
        "Fibro-reg NPCs",
        "Hom/Reg-like NPCs"
      )
    )
  )

readr::write_csv(
  score_df,
  file.path(table_dir, "clean_npc_ferroptosis_driver_suppressor_scores_per_cell.csv")
)

gene_set_summary <- tibble(
  score = c("FDS", "FSS"),
  ferrdb_category = c("ferroptosis driver", "ferroptosis suppressor"),
  ferrdb_genes_total_human_filtered = c(length(driver_genes_all), length(suppressor_genes_all)),
  genes_found_in_expression_matrix = c(length(driver_genes), length(suppressor_genes)),
  genes_used = c(paste(driver_genes, collapse = ";"), paste(suppressor_genes, collapse = ";"))
)
readr::write_csv(gene_set_summary, file.path(table_dir, "clean_npc_ferroptosis_score_gene_sets.csv"))

score_summary <- score_df %>%
  pivot_longer(c(FDS, FSS), names_to = "score", values_to = "value") %>%
  group_by(score, npc_subtype_final, article_group) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(value, na.rm = TRUE),
    median_score = median(value, na.rm = TRUE),
    sd_score = sd(value, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(score_summary, file.path(table_dir, "clean_npc_ferroptosis_score_summary_by_subtype_group.csv"))

score_stats <- score_df %>%
  pivot_longer(c(FDS, FSS), names_to = "score", values_to = "value") %>%
  group_by(score, npc_subtype_final) %>%
  summarise(
    n_mild = sum(article_group == "mild" & !is.na(value)),
    n_severe = sum(article_group == "severe" & !is.na(value)),
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
readr::write_csv(score_stats, file.path(table_dir, "clean_npc_ferroptosis_score_wilcox_stats.csv"))

plot_score_violin <- function(df, score_col, y_label, title, output_file) {
  y_limits <- df %>%
    group_by(npc_subtype_final) %>%
    summarise(y = max(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y + 0.08 * diff(range(df[[score_col]], na.rm = TRUE)))

  labels <- score_stats %>%
    filter(score == score_col) %>%
    select(npc_subtype_final, p_label) %>%
    left_join(y_limits, by = "npc_subtype_final")

  ggplot(df, aes(x = npc_subtype_final, y = .data[[score_col]], fill = article_group)) +
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
      data = labels,
      aes(x = npc_subtype_final, y = y, label = p_label),
      inherit.aes = FALSE,
      size = 4.2
    ) +
    scale_fill_manual(values = c(mild = "#4E79A7", severe = "#E15759")) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      fill = "Degeneration group"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "top"
    )

  ggsave(output_file, width = 10.5, height = 5.8, dpi = 300)
}

p_fds <- plot_score_violin(
  score_df,
  "FDS",
  "Ferroptosis driver score",
  "Fig. 2c-style: FerrDb driver score across NPC subtypes",
  file.path(figure_dir, "25_clean_npc_ferroptosis_driver_score_violin.png")
)

p_fss <- plot_score_violin(
  score_df,
  "FSS",
  "Ferroptosis suppressor score",
  "Fig. 2d-style: FerrDb suppressor score across NPC subtypes",
  file.path(figure_dir, "26_clean_npc_ferroptosis_suppressor_score_violin.png")
)

score_long <- score_df %>%
  pivot_longer(c(FDS, FSS), names_to = "score", values_to = "value")

p_compare <- ggplot(score_long, aes(x = article_group, y = value, fill = article_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.86, linewidth = 0.25) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.72, linewidth = 0.25) +
  facet_grid(score ~ npc_subtype_final, scales = "free_y") +
  scale_fill_manual(values = c(mild = "#4E79A7", severe = "#E15759")) +
  labs(
    title = "FerrDb ferroptosis scores: mild vs severe",
    x = NULL,
    y = "AddModuleScore value",
    fill = "Degeneration group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
ggsave(file.path(figure_dir, "27_clean_npc_ferroptosis_scores_mild_vs_severe_facets.png"), p_compare, width = 13, height = 6.8, dpi = 300)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("FerrDb FDS/FSS score plots finished.")
