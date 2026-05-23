suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(patchwork)
  library(Matrix)
})

project_root <- normalizePath(".", winslash = "/")
raw_dir <- file.path(project_root, "data/raw/GSE165722_RAW")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures")
table_dir <- file.path(project_root, "results/tables")
tmp_dir <- file.path(project_root, "tmp")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)
options(future.globals.maxSize = 8 * 1024^3)

set.seed(20260517)

metadata_path <- file.path(project_root, "data/raw/GSE165722_sample_metadata.csv")
sample_metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)

message("Listing raw count files...")
count_files <- list.files(raw_dir, pattern = "counts\\.tsv\\.gz$", full.names = TRUE)
print(basename(count_files))

read_sample_counts <- function(count_file) {
  sample_prefix <- sub("\\.counts\\.tsv\\.gz$", "", basename(count_file))
  sample_id <- sub("^GSM[0-9]+_", "", sample_prefix)
  cellname_file <- file.path(raw_dir, paste0(sample_prefix, ".cellname.txt.gz"))

  message("Reading normalized count table: ", basename(count_file))
  counts_tbl <- readr::read_tsv(count_file, show_col_types = FALSE)
  gene_col <- names(counts_tbl)[1]
  mat <- counts_tbl %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()
  storage.mode(mat) <- "numeric"

  if (file.exists(cellname_file)) {
    cell_names_tbl <- readr::read_tsv(cellname_file, show_col_types = FALSE)
    if ("CellName" %in% colnames(cell_names_tbl) && nrow(cell_names_tbl) == ncol(mat)) {
      colnames(mat) <- paste(sample_id, cell_names_tbl$CellName, sep = "_")
    } else if ("CellIndex" %in% colnames(cell_names_tbl) && nrow(cell_names_tbl) == ncol(mat)) {
      colnames(mat) <- paste(sample_id, cell_names_tbl$CellIndex, sep = "_")
    } else {
      warning("Cell-name count does not match matrix columns for ", sample_id, ". Keeping original column names.")
      colnames(mat) <- paste(sample_id, colnames(mat), sep = "_")
    }
  } else {
    colnames(mat) <- paste(sample_id, colnames(mat), sep = "_")
  }

  Matrix::Matrix(mat, sparse = TRUE)
}

message("Creating Seurat objects...")
objects <- lapply(count_files, function(path) {
  sample_prefix <- sub("\\.counts\\.tsv\\.gz$", "", basename(path))
  geo_accession <- sub("_.*$", "", sample_prefix)
  sample_id <- sub("^GSM[0-9]+_", "", sample_prefix)
  sample_info <- sample_metadata %>% filter(.data$geo_accession == !!geo_accession)
  counts <- read_sample_counts(path)
  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE165722",
    min.cells = 3,
    min.features = 200
  )
  obj$sample_id <- sample_id
  obj$geo_accession <- geo_accession
  obj$grade <- sample_info$grade[1]
  obj$grade_num <- sample_info$grade_num[1]
  obj$degeneration_group <- sample_info$degeneration_group[1]
  obj
})

names(objects) <- vapply(objects, function(x) unique(x$sample_id), character(1))

message("Merging samples...")
idd <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "GSE165722")
idd[["percent.mt"]] <- PercentageFeatureSet(idd, pattern = "^MT-")

qc_plot <- VlnPlot(
  idd,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  pt.size = 0,
  ncol = 3
)
ggsave(file.path(figure_dir, "01_qc_violin_by_sample.png"), qc_plot, width = 12, height = 5, dpi = 300)

message("Filtering cells with paper-like QC thresholds...")
idd <- subset(
  idd,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < median(nFeature_RNA) * 2 &
    percent.mt < 20
)

message("Normalizing, PCA, Harmony integration, clustering, UMAP...")
idd <- NormalizeData(idd)
idd <- FindVariableFeatures(idd, selection.method = "vst", nfeatures = 2000)
idd <- ScaleData(idd, features = rownames(idd))
idd <- RunPCA(idd, features = VariableFeatures(idd), npcs = 30)
idd <- RunHarmony(
  object = idd,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:30
)
idd <- FindNeighbors(idd, reduction = "harmony", dims = 1:30)
idd <- FindClusters(idd, resolution = 0.4)
idd <- RunUMAP(idd, reduction = "harmony", dims = 1:30)

saveRDS(idd, file.path(processed_dir, "GSE165722_seurat_harmony_clustered_preannotation.rds"))

p_umap_cluster <- DimPlot(idd, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("GSE165722 clusters")
p_umap_sample <- DimPlot(idd, reduction = "umap", group.by = "sample_id") +
  ggtitle("GSE165722 samples")
p_umap_grade <- DimPlot(idd, reduction = "umap", group.by = "grade") +
  ggtitle("GSE165722 degeneration grades")

ggsave(file.path(figure_dir, "02_umap_clusters.png"), p_umap_cluster, width = 7, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "03_umap_samples.png"), p_umap_sample, width = 7, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "04_umap_grades.png"), p_umap_grade, width = 7, height = 6, dpi = 300)

message("Finding cluster markers...")
markers <- FindAllMarkers(idd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
readr::write_csv(markers, file.path(table_dir, "cluster_markers.csv"))

if (nrow(markers) > 0 && "cluster" %in% colnames(markers)) {
  fc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(markers) ~ "avg_log2FC",
    "avg_logFC" %in% colnames(markers) ~ "avg_logFC",
    TRUE ~ NA_character_
  )

  if (!is.na(fc_col)) {
    top_markers <- markers %>%
      group_by(.data$cluster) %>%
      slice_max(.data[[fc_col]], n = 10, with_ties = FALSE) %>%
      ungroup()
    readr::write_csv(top_markers, file.path(table_dir, "top10_cluster_markers.csv"))
  } else {
    warning("Marker table has no avg_log2FC/avg_logFC column. Skipping top marker table.")
  }
} else {
  warning("FindAllMarkers returned no marker rows. Skipping top marker table.")
}

major_markers <- c(
  "COL2A1", "ACAN", "SOX9", "KRT19", "MMP3", "MMP13",
  "COL1A1", "COL3A1", "PTPRC", "LYZ", "PECAM1", "VWF"
)
major_markers <- major_markers[major_markers %in% rownames(idd)]

if (length(major_markers) > 0) {
  p_features <- FeaturePlot(idd, features = major_markers, ncol = 3)
  ggsave(file.path(figure_dir, "05_major_marker_featureplots.png"), p_features, width = 12, height = 10, dpi = 300)
  p_dot <- DotPlot(idd, features = major_markers, group.by = "seurat_clusters") +
    RotatedAxis() +
    ggtitle("Major marker genes by cluster")
  ggsave(file.path(figure_dir, "05b_major_marker_dotplot.png"), p_dot, width = 11, height = 6, dpi = 300)
}

message("Adding exploratory gene-set scores...")
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
  p_scores <- VlnPlot(idd, features = score_cols, group.by = "seurat_clusters", pt.size = 0, ncol = 1)
  p_scores_grade <- VlnPlot(idd, features = score_cols, group.by = "grade", pt.size = 0, ncol = 1)
  ggsave(file.path(figure_dir, "06_gene_set_scores_by_cluster.png"), p_scores, width = 9, height = 10, dpi = 300)
  ggsave(file.path(figure_dir, "07_gene_set_scores_by_grade.png"), p_scores_grade, width = 9, height = 10, dpi = 300)
}

saveRDS(idd, file.path(processed_dir, "GSE165722_seurat_harmony_clustered.rds"))

unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Finished. Outputs are in results/figures, results/tables, and data/processed.")
