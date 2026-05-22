# R version 4.4.3
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("Seurat")  # 5.3.0
BiocManager::install(c(
  "SoupX",           # 1.6.2
  "scDblFinder",     # 1.18.0
  "SingleCellExperiment", # 1.26.0
  "BiocParallel",    # 1.38.0
  "scran",           # 1.32.0
  "scry"             # 1.16.0
))
