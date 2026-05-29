# Install all R dependencies for pseudotime_02_tradeseq.R
# Run once:  Rscript pseudotime_00_install_r_deps.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "zellkonverter",         # read h5ad in R
  "SingleCellExperiment",  # data container
  "slingshot",             # trajectory inference
  "tradeSeq",              # DE along pseudotime
  "condiments",            # multi-condition trajectory comparison
  "BiocParallel"          # parallel fitting
)

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "viridis", "patchwork", "ggrepel")

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

cat("All packages installed.\n")
