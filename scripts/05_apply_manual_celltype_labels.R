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

clustered_rds <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered.rds")
labeled_rds <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered_labeled_v1.rds")
idd <- if (file.exists(clustered_rds)) {
  readRDS(clustered_rds)
} else {
  readRDS(labeled_rds)
}

cluster_labels <- c(
  "0" = "Activated myeloid/neutrophil",
  "1" = "NP/chondrocyte-like 1",
  "2" = "Inflammatory myeloid",
  "3" = "Fibro-ECM/stromal",
  "4" = "NP/chondrocyte-like 2",
  "5" = "Neutrophil granule-high",
  "6" = "T/NK cells",
  "7" = "Inflammatory neutrophil",
  "8" = "Immature neutrophil",
  "9" = "Inflammatory monocyte/macrophage",
  "10" = "Erythroid cells",
  "11" = "C1Q macrophages",
  "12" = "B/plasma cells",
  "13" = "Erythroid/platelet-like"
)

idd$manual_celltype_v1 <- unname(cluster_labels[as.character(idd$seurat_clusters)])

annotation_table <- tibble(
  cluster = names(cluster_labels),
  manual_celltype_v1 = unname(cluster_labels),
  rationale = c(
    "CXCR2, TLR2, ITGAX, PPIF; inflammatory myeloid/neutrophil-like signal",
    "FRZB, WIF1, S100B, CLEC3A, SCRG1; NP/chondrocyte-like homeostatic markers",
    "IL1RN, CCL20, VNN2, IL18R1; inflammatory myeloid signal",
    "POSTN, COL1A1, COL3A1, MMP2, TNC; collagen/ECM stromal phenotype",
    "COL2A1, ACAN, COL11A1, PRG4, VCAN; NP/chondrocyte-like ECM phenotype",
    "CAMP, LTF, MMP8, LCN2, OLFM4; mature neutrophil granule genes",
    "CD8A, GNLY, GZMA, CD3D, CD3G; T/NK lineage markers",
    "MMP9, S100A12, CD177, PADI4; inflammatory neutrophil markers",
    "DEFA4, CTSG, MPO, ELANE, MS4A3; immature granulocyte/neutrophil markers",
    "EREG, CXCL3, HLA-DRA, IL1A, CCL3; activated inflammatory monocyte/macrophage signal",
    "CA1, CA2, AHSP, HEMGN, GYPA; erythroid lineage markers",
    "C1QA, C1QB, C1QC, MSR1, APOC1; macrophage/complement markers",
    "JCHAIN, IGKC, IGHM, CD79A, MS4A1; B cell/plasma cell markers",
    "TRIM58, SLFN14, TSPO2, GCNT2; erythroid/platelet-like low-frequency cluster"
  )
)
readr::write_csv(annotation_table, file.path(table_dir, "manual_celltype_annotation_v1.csv"))

p_celltype <- DimPlot(
  idd,
  reduction = "umap",
  group.by = "manual_celltype_v1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("GSE165722 cell-type annotation")

ggsave(file.path(figure_dir, "08_umap_manual_celltypes_v1.png"), p_celltype, width = 10, height = 7, dpi = 300)

npc_candidate_clusters <- c("1", "3", "4")
idd$npc_candidate_v1 <- ifelse(
  as.character(idd$seurat_clusters) %in% npc_candidate_clusters,
  "NPC/stromal candidate",
  "non-NPC candidate"
)

p_npc_candidate <- DimPlot(
  idd,
  reduction = "umap",
  group.by = "npc_candidate_v1"
) +
  ggtitle("NPC/stromal candidate clusters for next-step subclustering")

ggsave(file.path(figure_dir, "09_umap_npc_candidate_v1.png"), p_npc_candidate, width = 8, height = 6, dpi = 300)

saveRDS(idd, file.path(processed_dir, "GSE165722_seurat_harmony_clustered_labeled_v1.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Manual cell type labels v1 applied.")
