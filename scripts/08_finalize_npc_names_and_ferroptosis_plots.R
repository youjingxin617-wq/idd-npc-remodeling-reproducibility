suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures")
table_dir <- file.path(project_root, "results/tables")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

clean_npc_v1 <- file.path(processed_dir, "GSE165722_clean_npc_subset_v1.rds")
clean_npc_final <- file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds")
clean_npc <- if (file.exists(clean_npc_v1)) {
  readRDS(clean_npc_v1)
} else {
  readRDS(clean_npc_final)
}
DefaultAssay(clean_npc) <- "RNA"

final_labels <- c(
  "0" = "Hom/Reg-like NPCs",
  "1" = "ECM/Adh-NPCs",
  "2" = "Hom-NPCs",
  "3" = "Fibro-NPCs",
  "4" = "Fibro-reg NPCs",
  "5" = "Eff-NPCs",
  "7" = "Ht-NPCs"
)

clean_npc$npc_subtype_final <- factor(
  unname(final_labels[as.character(clean_npc$npc_clean_cluster)]),
  levels = c("ECM/Adh-NPCs", "Hom/Reg-like NPCs", "Hom-NPCs", "Eff-NPCs", "Ht-NPCs", "Fibro-NPCs", "Fibro-reg NPCs")
)

final_annotation <- tibble(
  npc_cluster = names(final_labels),
  npc_subtype_final = unname(final_labels),
  key_markers = c(
    "CHRDL2, FRZB, WIF1, IGF2, CLEC3A",
    "PRG4, COL2A1, ACAN, COL11A1, VCAN, FN1",
    "SCRG1, MGP, C2orf40, MT1G",
    "POSTN, COL1A1, COL3A1, TNC, COL6A3",
    "COL12A1, MMP2, CRTAC1, SOCS3, CDH11",
    "CCL20, CXCL2, CHI3L1, CHI3L2, MMP3, SOD2",
    "IBSP, COL10A1, CYTL1, COL2A1, COL11A1"
  ),
  naming_logic = c(
    "Homeostatic/regulatory matrix signals supported by CHRDL2, FRZB, WIF1, IGF2, and CLEC3A",
    "ECM structural and adhesion-associated genes; closest to ECM-reg/Adh-NPC state",
    "Homeostatic NPC markers and low inflammatory/fibrotic signal",
    "Strong collagen/fibrotic ECM remodeling signature",
    "Fibro-ECM plus regulatory/stress-response markers",
    "Inflammatory and catabolic effector genes, including MMP3 and cytokine-response genes",
    "Hypertrophic/chondrocyte-like markers such as IBSP and COL10A1"
  )
)
readr::write_csv(final_annotation, file.path(table_dir, "npc_final_subtype_annotation_v2.csv"))

p_final_tsne <- DimPlot(
  clean_npc,
  reduction = "tsne",
  group.by = "npc_subtype_final",
  split.by = "article_group"
) +
  ggtitle("NPC subtype composition in mild and severe degeneration")
ggsave(file.path(figure_dir, "20_final_npc_tsne_subtypes_mild_vs_severe.png"), p_final_tsne, width = 12, height = 5.8, dpi = 300)

ferroptosis_genes <- c("ACSL4", "SLC39A14", "CHAC1", "GPX4", "FTH1", "HSPB1")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(clean_npc))

if (length(ferroptosis_genes) > 0) {
  p_vln <- VlnPlot(
    clean_npc,
    features = ferroptosis_genes,
    group.by = "npc_subtype_final",
    split.by = "article_group",
    pt.size = 0,
    ncol = 2,
    same.y.lims = FALSE
  ) &
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top"
    )
  ggsave(file.path(figure_dir, "21_clean_npc_ferroptosis_gene_violins_named.png"), p_vln, width = 13, height = 10, dpi = 300)

  p_dot_subtype <- DotPlot(clean_npc, features = ferroptosis_genes, group.by = "npc_subtype_final") +
    RotatedAxis() +
    ggtitle("Ferroptosis-related genes across final NPC subtypes")
  ggsave(file.path(figure_dir, "22_clean_npc_ferroptosis_gene_dotplot_by_subtype.png"), p_dot_subtype, width = 9, height = 5, dpi = 300)

  clean_npc$subtype_group <- paste(clean_npc$npc_subtype_final, clean_npc$article_group, sep = "_")
  p_dot_group <- DotPlot(clean_npc, features = ferroptosis_genes, group.by = "subtype_group") +
    RotatedAxis() +
    ggtitle("Ferroptosis-related genes by NPC subtype and degeneration group")
  ggsave(file.path(figure_dir, "23_clean_npc_ferroptosis_gene_dotplot_by_subtype_group.png"), p_dot_group, width = 13, height = 6, dpi = 300)

  avg_expr <- AverageExpression(
    clean_npc,
    features = ferroptosis_genes,
    group.by = "subtype_group",
    assays = "RNA",
    slot = "data"
  )$RNA
  readr::write_csv(
    as.data.frame(avg_expr) %>% rownames_to_column("gene"),
    file.path(table_dir, "clean_npc_ferroptosis_average_expression.csv")
  )

  png(file.path(figure_dir, "24_clean_npc_ferroptosis_average_expression_heatmap.png"), width = 1800, height = 900, res = 180)
  pheatmap(
    avg_expr,
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize = 9,
    main = "Average ferroptosis-related gene expression"
  )
  dev.off()
}

readr::write_csv(
  clean_npc@meta.data %>%
    count(npc_subtype_final, article_group, name = "n_cells") %>%
    group_by(article_group) %>%
    mutate(proportion = n_cells / sum(n_cells)) %>%
    ungroup(),
  file.path(table_dir, "npc_final_subtype_counts_and_proportions.csv")
)

saveRDS(clean_npc, file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Final NPC labels and optimized ferroptosis plots finished.")
