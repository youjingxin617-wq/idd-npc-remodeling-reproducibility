options(timeout = 1200)
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  download.file.method = "libcurl"
)

cran_packages <- c(
  "Seurat",
  "harmony",
  "tidyverse",
  "patchwork",
  "Matrix",
  "pheatmap",
  "msigdbr"
)

bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(cran_packages, install_if_missing))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::repositories()

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

missing_packages <- c(cran_packages, bioc_packages)[
  !vapply(c(cran_packages, bioc_packages), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  message("These packages are still missing: ", paste(missing_packages, collapse = ", "))
  message("Run scripts/00b_install_failed_packages.R to retry the slow Bioconductor packages.")
} else {
  message("All required packages are installed.")
}

message("Package installation check finished.")
