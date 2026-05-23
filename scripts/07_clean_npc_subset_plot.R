suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
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

npc <- readRDS(file.path(processed_dir, "GSE165722_npc_subset_reclustered.rds"))
Idents(npc) <- "npc_cluster"

# Remove obvious non-NPC contaminants after NPC/stromal candidate reclustering:
# cluster 6: LYZ/S100A8/IL1B myeloid
# cluster 8: PECAM1/KDR/EMCN endothelial
clean_npc <- subset(npc, idents = c("6", "8"), invert = TRUE)
clean_npc$npc_clean_cluster <- droplevels(clean_npc$npc_cluster)

clean_labels <- c(
  "0" = "Hom/Reg-like NPC",
  "1" = "ECM/Adh-like NPC",
  "2" = "Homeostatic NPC",
  "3" = "Fibro-NPC",
  "4" = "Fibro/Reg-like NPC",
  "5" = "Inflammatory Eff-NPC",
  "7" = "Hypertrophic-like NPC"
)
clean_npc$npc_subtype_v1 <- unname(clean_labels[as.character(clean_npc$npc_clean_cluster)])

label_table <- tibble(
  npc_cluster = names(clean_labels),
  npc_subtype_v1 = unname(clean_labels),
  rationale = c(
    "CHRDL2, FRZB, WIF1, IGF2, CLEC3A suggest homeostatic/regulatory NPC state",
    "PRG4, COL2A1, ACAN, COL11A1, VCAN, FN1 suggest ECM/adhesion NPC state",
    "SCRG1, MGP, C2orf40 and mitochondrial/ribosomal genes suggest homeostatic NPC state",
    "POSTN, COL1A1, COL3A1, TNC, COL6A3 suggest fibrotic ECM-rich NPC state",
    "COL12A1, MMP2, CRTAC1, SOCS3, CDH11 suggest fibro/regulatory ECM state",
    "CCL20, CXCL2, CHI3L1, CHI3L2, MMP3, SOD2 suggest inflammatory effector NPC state",
    "IBSP, COL10A1, CYTL1, COL2A1, COL11A1 suggest hypertrophic/chondrocyte-like NPC state"
  )
)
readr::write_csv(label_table, file.path(table_dir, "npc_clean_subtype_annotation_v1.csv"))

p_tsne_named <- DimPlot(
  clean_npc,
  reduction = "tsne",
  group.by = "npc_subtype_v1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Clean NPC-like subtypes v1")

p_tsne_named_split <- DimPlot(
  clean_npc,
  reduction = "tsne",
  group.by = "npc_subtype_v1",
  split.by = "article_group"
) +
  ggtitle("Clean NPC-like subtypes: mild vs severe")

p_cluster_split <- DimPlot(
  clean_npc,
  reduction = "tsne",
  group.by = "npc_clean_cluster",
  split.by = "article_group",
  label = TRUE
) +
  ggtitle("Clean NPC-like clusters: mild vs severe")

ggsave(file.path(figure_dir, "16_clean_npc_tsne_subtype_names.png"), p_tsne_named, width = 10, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "17_clean_npc_tsne_subtype_names_mild_vs_severe.png"), p_tsne_named_split, width = 12, height = 5.8, dpi = 300)
ggsave(file.path(figure_dir, "18_clean_npc_tsne_clusters_mild_vs_severe.png"), p_cluster_split, width = 12, height = 5.8, dpi = 300)

article_marker_genes <- c(
  "NFKBIA", "NFKBIZ", "MMP3", "CHI3L2", "CHI3L1", "TNC",
  "COL1A2", "COL3A1", "COL1A1", "COL14A1", "COL12A1",
  "MMP2", "MSMO1", "PRDX4", "ITM2A", "FGFBP2", "SCRG1",
  "AOC2", "COL15A1", "CRTAC1", "FMOD", "FN1", "KLF4",
  "SOX9", "ACAN", "COL2A1", "IGF2", "CHRDL2", "CILP2",
  "SEMA3A", "SLPI", "VCAN"
)
article_marker_genes <- intersect(article_marker_genes, rownames(clean_npc))

if (length(article_marker_genes) > 0) {
  p_dot <- DotPlot(clean_npc, features = article_marker_genes, group.by = "npc_subtype_v1") +
    RotatedAxis() +
    ggtitle("Article marker genes across clean NPC-like subtypes")
  ggsave(file.path(figure_dir, "19_clean_npc_article_marker_dotplot_named.png"), p_dot, width = 13, height = 8, dpi = 300)
}

readr::write_csv(
  clean_npc@meta.data %>%
    count(npc_subtype_v1, article_group, name = "n_cells"),
  file.path(table_dir, "clean_npc_subtype_counts_mild_severe.csv")
)

saveRDS(clean_npc, file.path(processed_dir, "GSE165722_clean_npc_subset_v1.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Clean NPC subset plots finished.")
