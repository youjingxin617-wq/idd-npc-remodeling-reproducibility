suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
})

project_root <- normalizePath(".", winslash = "/")
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

set.seed(20260519)

message("Loading labeled whole-cell object...")
whole_path <- file.path(processed_dir, "GSE165722_seurat_harmony_clustered_labeled_v1.rds")
idd <- readRDS(whole_path)
DefaultAssay(idd) <- "RNA"

message("Joining Seurat v5 layers for subset re-analysis...")
idd <- JoinLayers(idd)

# Article-like setup:
# The paper analyzes human NP tissues with mild degeneration (grades II/III)
# versus severe degeneration (grades IV/V). Public GSE165722 metadata here has
# Grade I-IV samples, so Grade I is excluded and Grade IV is used as severe.
npc_candidate_clusters <- c("1", "3", "4")
Idents(idd) <- "seurat_clusters"
npc <- subset(
  idd,
  idents = npc_candidate_clusters,
  subset = grade_num >= 2
)

npc$article_group <- case_when(
  npc$grade_num %in% c(2, 3) ~ "mild",
  npc$grade_num >= 4 ~ "severe",
  TRUE ~ NA_character_
)
npc$article_group <- factor(npc$article_group, levels = c("mild", "severe"))
npc$source_whole_cluster <- as.character(npc$seurat_clusters)

readr::write_csv(
  npc@meta.data %>%
    count(sample_id, grade, grade_num, article_group, source_whole_cluster, name = "n_cells"),
  file.path(table_dir, "npc_subset_cell_counts_by_sample_cluster.csv")
)

message("Reprocessing NPC/stromal candidate subset...")
npc <- NormalizeData(npc, verbose = FALSE)
npc <- FindVariableFeatures(npc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
npc <- ScaleData(npc, features = rownames(npc), verbose = FALSE)
npc <- RunPCA(npc, features = VariableFeatures(npc), npcs = 30, verbose = FALSE)
npc <- RunHarmony(
  object = npc,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = 1:30,
  verbose = FALSE
)
npc <- FindNeighbors(npc, reduction = "harmony", dims = 1:30, verbose = FALSE)
npc <- FindClusters(npc, resolution = 0.4, verbose = FALSE)
npc <- RunUMAP(npc, reduction = "harmony", dims = 1:30, verbose = FALSE)
npc <- RunTSNE(npc, reduction = "harmony", dims = 1:30, check_duplicates = FALSE)

npc$npc_cluster <- factor(npc$seurat_clusters)

message("Drawing article-like NPC tSNE/UMAP panels...")
p_tsne_clusters <- DimPlot(npc, reduction = "tsne", group.by = "npc_cluster", label = TRUE, repel = TRUE) +
  ggtitle("NPC/stromal candidate subclusters")
p_tsne_split <- DimPlot(npc, reduction = "tsne", group.by = "npc_cluster", split.by = "article_group", label = TRUE) +
  ggtitle("NPC/stromal candidate subclusters: mild vs severe")
p_tsne_group <- DimPlot(npc, reduction = "tsne", group.by = "article_group") +
  ggtitle("NPC/stromal candidate cells by article-like group")
p_umap_split <- DimPlot(npc, reduction = "umap", group.by = "npc_cluster", split.by = "article_group", label = TRUE) +
  ggtitle("NPC/stromal candidate UMAP: mild vs severe")

ggsave(file.path(figure_dir, "10_npc_subset_tsne_clusters.png"), p_tsne_clusters, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "11_npc_subset_tsne_mild_vs_severe.png"), p_tsne_split, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(figure_dir, "12_npc_subset_tsne_article_group.png"), p_tsne_group, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "13_npc_subset_umap_mild_vs_severe.png"), p_umap_split, width = 11, height = 5.5, dpi = 300)

article_marker_genes <- c(
  "NFKBIA", "NFKBIZ", "MMP3", "CHI3L2", "CHI3L1", "TNC",
  "COL1A2", "COL3A1", "COL1A1", "COL14A1", "COL12A1",
  "MMP2", "MSMO1", "PRDX4", "ITM2A", "FGFBP2", "SCRG1",
  "AOC2", "COL15A1", "CRTAC1", "FMOD", "FN1", "KLF4",
  "SOX9", "ACAN", "COL2A1", "IGF2", "CHRDL2", "CILP2",
  "SEMA3A", "SLPI", "VCAN"
)
article_marker_genes <- intersect(article_marker_genes, rownames(npc))

if (length(article_marker_genes) > 0) {
  p_dot <- DotPlot(npc, features = article_marker_genes, group.by = "npc_cluster") +
    RotatedAxis() +
    ggtitle("Article marker genes across NPC/stromal candidate subclusters")
  ggsave(file.path(figure_dir, "14_npc_subset_article_marker_dotplot.png"), p_dot, width = 12, height = 8, dpi = 300)
}

ferroptosis_genes <- c("ACSL4", "SLC39A14", "CHAC1", "GPX4", "FTH1", "HSPB1")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(npc))

if (length(ferroptosis_genes) > 0) {
  p_ferro <- VlnPlot(
    npc,
    features = ferroptosis_genes,
    group.by = "npc_cluster",
    split.by = "article_group",
    pt.size = 0,
    ncol = 2
  )
  ggsave(file.path(figure_dir, "15_npc_subset_ferroptosis_gene_violins.png"), p_ferro, width = 12, height = 10, dpi = 300)
}

message("Finding NPC subset markers...")
Idents(npc) <- "npc_cluster"
markers <- FindAllMarkers(npc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
readr::write_csv(markers, file.path(table_dir, "npc_subset_cluster_markers.csv"))

if (nrow(markers) > 0) {
  fc_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
  top_markers <- markers %>%
    group_by(.data$cluster) %>%
    slice_max(.data[[fc_col]], n = 15, with_ties = FALSE) %>%
    ungroup()
  readr::write_csv(top_markers, file.path(table_dir, "npc_subset_top15_cluster_markers.csv"))
}

readr::write_csv(
  npc@meta.data %>%
    count(npc_cluster, article_group, source_whole_cluster, name = "n_cells"),
  file.path(table_dir, "npc_subset_cluster_composition_mild_severe.csv")
)

saveRDS(npc, file.path(processed_dir, "GSE165722_npc_subset_reclustered.rds"))
unlink(tmp_dir, recursive = TRUE, force = TRUE)

message("NPC subset reproduction finished.")
