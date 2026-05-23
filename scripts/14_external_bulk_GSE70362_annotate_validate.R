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
  "COL1A1", "COL1A2", "COL3A1", "COL6A2", "COMP", "BGN", "SPARC", "LRP1", "FN1", "POSTN",
  "GPX3", "CLU", "SLC39A14", "GPX4", "SLC7A11", "FTH1", "FTL", "HMOX1", "ACSL4",
  "ACAN", "COL2A1", "MMP3", "MMP13", "ADAMTS5"
)

expr <- readRDS(file.path(processed_dir, "GSE70362_expr_raw_series_matrix.rds"))
pheno <- readr::read_csv(file.path(table_dir, "GSE70362_pheno_raw.csv"), show_col_types = FALSE)

gene_info_path <- file.path(raw_dir, "Homo_sapiens.gene_info.gz")
if (!file.exists(gene_info_path)) {
  stop("Missing NCBI gene_info file: ", gene_info_path)
}

gene_info <- readr::read_tsv(gzfile(gene_info_path), show_col_types = FALSE)
if ("#tax_id" %in% colnames(gene_info)) {
  gene_info <- gene_info %>% rename(tax_id = `#tax_id`)
}
gene_map <- gene_info %>%
  transmute(entrez_id = as.character(GeneID), symbol = Symbol, description = description) %>%
  distinct(entrez_id, .keep_all = TRUE)

probe_map <- tibble(
  probe_id = rownames(expr),
  entrez_id_raw = str_remove(rownames(expr), "_at$"),
  numeric_entrez = str_detect(entrez_id_raw, "^[0-9]+$")
) %>%
  left_join(gene_map, by = c("entrez_id_raw" = "entrez_id"))

readr::write_csv(probe_map, file.path(table_dir, "GSE70362_probe_to_symbol_mapping.csv"))

mapping_summary <- tibble(
  n_probes = nrow(probe_map),
  n_numeric_entrez_like = sum(probe_map$numeric_entrez),
  n_mapped_symbol = sum(!is.na(probe_map$symbol)),
  n_unmapped = sum(is.na(probe_map$symbol))
)
readr::write_csv(mapping_summary, file.path(table_dir, "GSE70362_probe_mapping_summary.csv"))
print(mapping_summary)

mapped <- probe_map %>% filter(!is.na(symbol), symbol != "")
expr_mapped <- expr[mapped$probe_id, , drop = FALSE]
rownames(expr_mapped) <- mapped$symbol

# Collapse duplicated symbols by mean expression. Brainarray should be near one-to-one, but we do not assume it.
expr_symbol <- rowsum(expr_mapped, group = rownames(expr_mapped), reorder = FALSE)
gene_counts <- table(rownames(expr_mapped))
expr_symbol <- sweep(expr_symbol, 1, as.numeric(gene_counts[rownames(expr_symbol)]), "/")

saveRDS(expr_symbol, file.path(processed_dir, "GSE70362_expr_gene_symbol_collapsed.rds"))

pheno2 <- pheno %>%
  mutate(
    tissue = characteristics_3_tissue,
    thompson_grade = characteristics_2_thompson.grade,
    degeneration_group = case_when(
      tissue == "Nucleus pulposus" & thompson_grade %in% c("I", "I-II", "II") ~ "low",
      tissue == "Nucleus pulposus" & thompson_grade %in% c("IV", "V") ~ "high",
      tissue == "Nucleus pulposus" & thompson_grade == "III" ~ "intermediate",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(pheno2, file.path(table_dir, "GSE70362_pheno_annotated.csv"))

np_pheno <- pheno2 %>%
  filter(tissue == "Nucleus pulposus") %>%
  mutate(use_for_low_high = degeneration_group %in% c("low", "high"))
readr::write_csv(np_pheno, file.path(table_dir, "GSE70362_NP_pheno_low_high_definition.csv"))

np_group_summary <- np_pheno %>% count(thompson_grade, degeneration_group, use_for_low_high)
readr::write_csv(np_group_summary, file.path(table_dir, "GSE70362_NP_group_summary.csv"))
print(np_group_summary)

available_core <- tibble(
  gene = core_genes,
  available_in_GSE70362 = gene %in% rownames(expr_symbol)
)
readr::write_csv(available_core, file.path(table_dir, "GSE70362_core_gene_availability.csv"))
print(available_core, n = Inf)

expr_np <- expr_symbol[, np_pheno$geo_accession, drop = FALSE]
saveRDS(expr_np, file.path(processed_dir, "GSE70362_NP_expr_gene_symbol_collapsed.rds"))

# Prepare a long table for available core genes in low/high NP samples.
usable_samples <- np_pheno %>% filter(use_for_low_high)
usable_genes <- available_core %>% filter(available_in_GSE70362) %>% pull(gene)
core_long <- expr_symbol[usable_genes, usable_samples$geo_accession, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "geo_accession", values_to = "expression") %>%
  left_join(usable_samples %>% select(geo_accession, title, thompson_grade, degeneration_group, batch = characteristics_4_batch), by = "geo_accession")
readr::write_csv(core_long, file.path(table_dir, "GSE70362_NP_low_high_core_gene_expression_long.csv"))

stats <- core_long %>%
  group_by(gene) %>%
  summarise(
    n_low = sum(degeneration_group == "low"),
    n_high = sum(degeneration_group == "high"),
    mean_low = mean(expression[degeneration_group == "low"], na.rm = TRUE),
    mean_high = mean(expression[degeneration_group == "high"], na.rm = TRUE),
    log2FC_high_vs_low = mean_high - mean_low,
    p_value = if_else(
      n_low >= 3 & n_high >= 3,
      wilcox.test(expression[degeneration_group == "high"], expression[degeneration_group == "low"], exact = FALSE)$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, p_value)
readr::write_csv(stats, file.path(table_dir, "GSE70362_NP_low_high_core_gene_stats.csv"))
print(stats, n = Inf)

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("GSE70362 annotation and core gene validation prep finished.")
