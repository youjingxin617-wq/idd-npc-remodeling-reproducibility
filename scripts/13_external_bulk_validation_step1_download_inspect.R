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

for (d in c(raw_dir, processed_dir, figure_dir, table_dir, tmp_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

series_matrix_url <- function(gse) {
  prefix <- substr(gse, 1, nchar(gse) - 3)
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", prefix, "nnn/", gse, "/matrix/", gse, "_series_matrix.txt.gz")
}

download_series_matrix <- function(gse) {
  dest <- file.path(raw_dir, paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(dest)) {
    message("Downloading ", gse, " series matrix...")
    download.file(series_matrix_url(gse), destfile = dest, mode = "wb", quiet = FALSE)
  } else {
    message("Using existing file: ", dest)
  }
  dest
}

parse_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_title <- strsplit(lines[str_starts(lines, "!Sample_title")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  sample_geo <- strsplit(lines[str_starts(lines, "!Sample_geo_accession")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  sample_source <- strsplit(lines[str_starts(lines, "!Sample_source_name_ch1")], "\t")[[1]][-1] %>% str_remove_all('^"|"$')
  characteristics_lines <- lines[str_starts(lines, "!Sample_characteristics_ch1")]

  pheno <- tibble(
    geo_accession = sample_geo,
    title = sample_title,
    source_name = sample_source
  )

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

inspect_gse <- function(gse) {
  path <- download_series_matrix(gse)
  parsed <- parse_series_matrix(path)
  expr <- parsed$expr
  pheno <- parsed$pheno

  readr::write_csv(pheno, file.path(table_dir, paste0(gse, "_pheno_raw.csv")))
  saveRDS(expr, file.path(processed_dir, paste0(gse, "_expr_raw_series_matrix.rds")))

  summary <- tibble(
    gse = gse,
    n_features = nrow(expr),
    n_samples = ncol(expr),
    min_expr = min(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    max_expr = max(expr, na.rm = TRUE),
    pheno_columns = paste(colnames(pheno), collapse = ";")
  )
  readr::write_csv(summary, file.path(table_dir, paste0(gse, "_series_matrix_summary.csv")))

  message("\n", gse, " summary:")
  print(summary)
  message("\n", gse, " phenotype preview:")
  print(pheno, n = min(12, nrow(pheno)), width = Inf)

  invisible(list(expr = expr, pheno = pheno, summary = summary))
}

res70362 <- inspect_gse("GSE70362")
res56081 <- inspect_gse("GSE56081")

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("External bulk download and inspection finished.")
