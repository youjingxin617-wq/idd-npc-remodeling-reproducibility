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

object_path <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered_preannotation.rds")
if (!file.exists(object_path)) {
  object_path <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered.rds")
}
if (!file.exists(object_path)) {
  stop("No saved Seurat object found. Run scripts/02_reproduce_GSE165722_scRNA.R first.")
}

idd <- readRDS(object_path)

major_markers <- c(
  "COL2A1", "ACAN", "SOX9", "KRT19", "MMP3", "MMP13",
  "COL1A1", "COL3A1", "PTPRC", "LYZ", "PECAM1", "VWF"
)
major_markers <- intersect(major_markers, rownames(idd))

if (length(major_markers) > 0) {
  p_features <- FeaturePlot(idd, features = major_markers, ncol = 3)
  ggsave(file.path(figure_dir, "05_major_marker_featureplots.png"), p_features, width = 12, height = 10, dpi = 300)

  p_dot <- DotPlot(idd, features = major_markers, group.by = "seurat_clusters") +
    RotatedAxis() +
    ggtitle("Major marker genes by cluster")
  ggsave(file.path(figure_dir, "05b_major_marker_dotplot.png"), p_dot, width = 11, height = 6, dpi = 300)
}

cluster_counts <- idd@meta.data %>%
  count(seurat_clusters, sample_id, grade, degeneration_group, name = "n_cells")
readr::write_csv(cluster_counts, file.path(table_dir, "cluster_cell_counts_by_sample_grade.csv"))

ros_genes <- c("SOD1", "SOD2", "CAT", "GPX1", "GPX3", "NQO1", "HMOX1", "TXN", "TXNRD1", "PRDX1", "PRDX2")
ferroptosis_genes <- c("GPX4", "SLC7A11", "ACSL4", "TFRC", "FTH1", "FTL", "NCOA4", "ALOX15", "HMOX1", "NFE2L2")
inflammation_genes <- c("IL1B", "IL6", "TNF", "CXCL8", "CCL2", "NFKB1", "NFKBIA", "PTGS2", "STAT3")

gene_sets <- list(
  ROS = intersect(ros_genes, rownames(idd)),
  Ferroptosis = intersect(ferroptosis_genes, rownames(idd)),
  Inflammation = intersect(inflammation_genes, rownames(idd))
)
gene_sets <- gene_sets[lengths(gene_sets) > 0]

if (length(gene_sets) > 0) {
  idd <- AddModuleScore(idd, features = gene_sets, name = names(gene_sets))
  score_cols <- grep("^(ROS|Ferroptosis|Inflammation)[0-9]+$", colnames(idd@meta.data), value = TRUE)
  p_scores_cluster <- VlnPlot(idd, features = score_cols, group.by = "seurat_clusters", pt.size = 0, ncol = 1)
  p_scores_grade <- VlnPlot(idd, features = score_cols, group.by = "grade", pt.size = 0, ncol = 1)
  ggsave(file.path(figure_dir, "06_gene_set_scores_by_cluster.png"), p_scores_cluster, width = 9, height = 10, dpi = 300)
  ggsave(file.path(figure_dir, "07_gene_set_scores_by_grade.png"), p_scores_grade, width = 9, height = 10, dpi = 300)
}

saveRDS(idd, file.path(processed_dir, "GSE165722_seurat_harmony_clustered_annot_ready.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Annotation support figures finished.")
