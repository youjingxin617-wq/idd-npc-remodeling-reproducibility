options(timeout = 1800)
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  download.file.method = "libcurl"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

retry_install <- function(pkg, tries = 3) {
  for (i in seq_len(tries)) {
    message("Installing ", pkg, " attempt ", i, " of ", tries)
    try(
      BiocManager::install(pkg, ask = FALSE, update = FALSE),
      silent = TRUE
    )
    if (requireNamespace(pkg, quietly = TRUE)) {
      message(pkg, " installed successfully.")
      return(invisible(TRUE))
    }
    Sys.sleep(5)
  }
  warning(pkg, " is still not installed.")
  invisible(FALSE)
}

retry_install("org.Hs.eg.db")
retry_install("clusterProfiler")

missing <- c("org.Hs.eg.db", "clusterProfiler")[
  !vapply(c("org.Hs.eg.db", "clusterProfiler"), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  message("Still missing: ", paste(missing, collapse = ", "))
  message("You can still run the first Seurat reproduction script. These packages are mainly needed later for GO/KEGG enrichment.")
} else {
  message("Bioconductor package retry finished successfully.")
}
