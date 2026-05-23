suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- normalizePath(".", winslash = "/")
bulk_dir <- file.path(project_root, "data/external_bulk")
raw_dir <- file.path(bulk_dir, "raw")
processed_dir <- file.path(bulk_dir, "processed")
figure_dir <- file.path(project_root, "results/figures_external_bulk")
table_dir <- file.path(project_root, "results/tables_external_bulk")
tmp_dir <- file.path(project_root, "tmp")

for (d in c(processed_dir, figure_dir, table_dir, tmp_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

core_genes <- c(
  "ACAN", "COL2A1", "COL6A2", "BGN", "GPX4", "SLC39A14", "FTH1", "FTL",
  "COL1A1", "COL1A2", "COL3A1", "POSTN", "ADAMTS5", "MMP13", "GPX3", "CLU",
  "COMP", "SPARC", "LRP1", "FN1", "SLC7A11", "HMOX1", "ACSL4"
)

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_title <- strsplit(lines[str_starts(lines, "!Sample_title")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  sample_geo <- strsplit(lines[str_starts(lines, "!Sample_geo_accession")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  sample_source <- strsplit(lines[str_starts(lines, "!Sample_source_name_ch1")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  characteristics_lines <- lines[str_starts(lines, "!Sample_characteristics_ch1")]

  pheno <- tibble(geo_accession = sample_geo, title = sample_title, source_name = sample_source)
  if (length(characteristics_lines) > 0) {
    for (i in seq_along(characteristics_lines)) {
      values <- strsplit(characteristics_lines[[i]], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
      key <- values[[1]] %>% str_split_fixed(":", 2) %>% .[, 1] %>% make.names()
      values_clean <- str_replace(values, "^[^:]+:\\s*", "")
      pheno[[paste0("characteristics_", i, "_", key)]] <- values_clean
    }
  }

  begin <- which(lines == "!series_matrix_table_begin") + 1
  end <- which(lines == "!series_matrix_table_end") - 1
  expr_lines <- lines[begin:end]
  expr <- readr::read_tsv(I(expr_lines), show_col_types = FALSE)
  expr <- as.data.frame(expr)
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  expr[] <- lapply(expr, as.numeric)
  list(expr = as.matrix(expr), pheno = pheno)
}

series_path <- file.path(raw_dir, "GSE147383_series_matrix.txt.gz")
annot_path <- file.path(raw_dir, "GSE147383_HG-U133_Plus_2.na28.annot.csv.gz")
if (!file.exists(series_path)) stop("Missing series matrix: ", series_path)
if (!file.exists(annot_path)) stop("Missing annotation file: ", annot_path)

parsed <- parse_series_matrix(series_path)
expr <- parsed$expr
pheno <- parsed$pheno

readr::write_csv(pheno, file.path(table_dir, "GSE147383_pheno_raw.csv"))
saveRDS(expr, file.path(processed_dir, "GSE147383_expr_raw_series_matrix.rds"))

message("GSE147383 phenotype:")
print(pheno, n = Inf, width = Inf)

# Affymetrix annotation files have comment lines before the real CSV header.
annot_lines <- readLines(gzfile(annot_path), warn = FALSE)
header_idx <- which(str_detect(annot_lines, '^"Probe Set ID"'))[1]
if (is.na(header_idx)) stop("Could not find Probe Set ID header in annotation file.")
annot <- readr::read_csv(I(annot_lines[header_idx:length(annot_lines)]), show_col_types = FALSE)

candidate_symbol_cols <- intersect(c("Gene Symbol", "gene_assignment", "Gene Assignment"), colnames(annot))
if (length(candidate_symbol_cols) == 0) {
  stop("No obvious gene symbol column found. Columns include: ", paste(colnames(annot), collapse = ";"))
}

symbol_col <- candidate_symbol_cols[1]
probe_map <- annot %>%
  transmute(
    probe_id = as.character(`Probe Set ID`),
    gene_symbol_raw = as.character(.data[[symbol_col]])
  ) %>%
  mutate(
    # Affymetrix annotation often separates multiple symbols by ///.
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1] %>% str_trim(),
    gene_symbol = na_if(gene_symbol, "---"),
    gene_symbol = na_if(gene_symbol, "")
  )

readr::write_csv(probe_map, file.path(table_dir, "GSE147383_probe_to_symbol_mapping.csv"))

mapping_summary <- tibble(
  n_expression_probes = nrow(expr),
  n_annotation_rows = nrow(probe_map),
  n_expr_probes_mapped = sum(rownames(expr) %in% probe_map$probe_id[!is.na(probe_map$gene_symbol)]),
  n_expr_probes_unmapped = nrow(expr) - n_expr_probes_mapped
)
readr::write_csv(mapping_summary, file.path(table_dir, "GSE147383_probe_mapping_summary.csv"))
print(mapping_summary)

mapped <- tibble(probe_id = rownames(expr)) %>%
  left_join(probe_map, by = "probe_id") %>%
  filter(!is.na(gene_symbol))
expr_mapped <- expr[mapped$probe_id, , drop = FALSE]
rownames(expr_mapped) <- mapped$gene_symbol
expr_symbol <- rowsum(expr_mapped, group = rownames(expr_mapped), reorder = FALSE)
gene_counts <- table(rownames(expr_mapped))
expr_symbol <- sweep(expr_symbol, 1, as.numeric(gene_counts[rownames(expr_symbol)]), "/")
saveRDS(expr_symbol, file.path(processed_dir, "GSE147383_expr_gene_symbol_collapsed.rds"))

pheno2 <- pheno %>%
  mutate(
    tissue = case_when(
      str_detect(str_to_lower(characteristics_4_tissue), "nucleus") | str_detect(str_to_lower(title), "np") | str_detect(str_to_lower(source_name), "nucleus") ~ "Nucleus pulposus",
      str_detect(str_to_lower(characteristics_4_tissue), "annulus") | str_detect(str_to_lower(title), "af") | str_detect(str_to_lower(source_name), "annulus") ~ "Annulus fibrosus",
      TRUE ~ source_name
    ),
    degeneration_group = case_when(
      str_detect(str_to_lower(characteristics_5_tissue.type), "non") | str_detect(str_to_lower(title), "young|non") | str_detect(str_to_lower(source_name), "young|non") ~ "young_non_degenerated",
      str_detect(str_to_lower(characteristics_5_tissue.type), "deg") | str_detect(str_to_lower(title), "aged|deg") | str_detect(str_to_lower(source_name), "aged|deg") ~ "aged_degenerated",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(pheno2, file.path(table_dir, "GSE147383_pheno_annotated.csv"))
print(pheno2, n = Inf, width = Inf)

available_core <- tibble(gene = core_genes, available_in_GSE147383 = gene %in% rownames(expr_symbol))
readr::write_csv(available_core, file.path(table_dir, "GSE147383_core_gene_availability.csv"))
print(available_core, n = Inf)

np_pheno <- pheno2 %>% filter(tissue == "Nucleus pulposus", !is.na(degeneration_group))
readr::write_csv(np_pheno, file.path(table_dir, "GSE147383_NP_pheno_group_definition.csv"))
print(np_pheno, n = Inf, width = Inf)

usable_genes <- available_core %>% filter(available_in_GSE147383) %>% pull(gene)
if (nrow(np_pheno) >= 4 && length(usable_genes) > 0) {
  long_gene <- expr_symbol[usable_genes, np_pheno$geo_accession, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "geo_accession", values_to = "expression") %>%
    left_join(np_pheno %>% select(geo_accession, title, tissue, degeneration_group), by = "geo_accession")
  readr::write_csv(long_gene, file.path(table_dir, "GSE147383_NP_core_gene_expression_long.csv"))

  stats <- long_gene %>%
    group_by(gene) %>%
    summarise(
      n_young = sum(degeneration_group == "young_non_degenerated"),
      n_aged = sum(degeneration_group == "aged_degenerated"),
      mean_young = mean(expression[degeneration_group == "young_non_degenerated"], na.rm = TRUE),
      mean_aged = mean(expression[degeneration_group == "aged_degenerated"], na.rm = TRUE),
      delta_aged_minus_young = mean_aged - mean_young,
      p_value = if_else(
        n_young >= 2 & n_aged >= 2,
        wilcox.test(expression[degeneration_group == "aged_degenerated"], expression[degeneration_group == "young_non_degenerated"], exact = FALSE)$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  readr::write_csv(stats, file.path(table_dir, "GSE147383_NP_core_gene_stats.csv"))
  print(stats, n = Inf)
} else {
  warning("NP pheno or usable genes insufficient for trend analysis.")
}

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("GSE147383 inspection finished.")
