suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
table_dir <- file.path(project_root, "results/tables")
tmp_dir <- file.path(project_root, "tmp")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

object_path <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered.rds")
idd <- readRDS(object_path)

# Seurat v5 stores each sample in separate assay layers after merge. JoinLayers
# is needed before marker testing across all cells.
idd <- JoinLayers(idd)
DefaultAssay(idd) <- "RNA"
Idents(idd) <- "seurat_clusters"

markers <- FindAllMarkers(
  idd,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

readr::write_csv(markers, file.path(table_dir, "cluster_markers.csv"))

fc_col <- dplyr::case_when(
  "avg_log2FC" %in% colnames(markers) ~ "avg_log2FC",
  "avg_logFC" %in% colnames(markers) ~ "avg_logFC",
  TRUE ~ NA_character_
)

if (nrow(markers) > 0 && !is.na(fc_col)) {
  top_markers <- markers %>%
    group_by(.data$cluster) %>%
    slice_max(.data[[fc_col]], n = 10, with_ties = FALSE) %>%
    ungroup()
  readr::write_csv(top_markers, file.path(table_dir, "top10_cluster_markers.csv"))

  annotation_hints <- top_markers %>%
    group_by(.data$cluster) %>%
    summarise(top_genes = paste(.data$gene, collapse = ", "), .groups = "drop") %>%
    mutate(
      annotation_guess = case_when(
        str_detect(top_genes, "ACAN|COL2A1|SOX9|KRT19|S100B|CLEC3A|ITM2A|CNMD") ~ "NP/chondrocyte-like cell candidate",
        str_detect(top_genes, "COL1A1|COL1A2|COL3A1|DCN|LUM|MMP2|FBN1") ~ "fibroblast/ECM cell candidate",
        str_detect(top_genes, "PTPRC|LYZ|LCP1|SRGN|S100A8|S100A9|CXCL8") ~ "immune/myeloid cell candidate",
        str_detect(top_genes, "NKG7|TRAC|TRBC|CD3D|CD8A|GZMA|GNLY") ~ "T/NK cell candidate",
        str_detect(top_genes, "PECAM1|VWF|KDR|CLDN5") ~ "endothelial cell candidate",
        str_detect(top_genes, "HBB|HBA1|HBA2|SLC4A1|ALAS2|GYPA") ~ "erythroid/RBC candidate",
        TRUE ~ "needs manual review"
      )
    )
  readr::write_csv(annotation_hints, file.path(table_dir, "cluster_annotation_hints.csv"))
}

saveRDS(idd, file.path(processed_dir, "GSE165722_seurat_harmony_clustered_joined.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("Marker detection finished.")
