packages <- read.csv(file.path('~/Documents/Research codes/cell_cell_communication/data/R_packages.csv'))[, -1]
base_packages <- as.data.frame(installed.packages()) 
to_install <- setdiff(packages$Package, base_packages$Package) 
install.packages(to_install)

# List of packages to install
packages <- c(
  "Matrix", "viridis", "harmony", "ggpubr", "tictoc", "RColorBrewer", "Hmisc",
  "corrplot", "grid", "gridExtra", "igraph", "ggrepel", "readxl", "conflicted",
  "dplyr", "parallel", "stringr", "Seurat", "SingleCellExperiment", "tidyr",
  "GSA", "limma", "tidyverse", "cowplot", "patchwork", "ggplot2", "WGCNA",
  "hdWGCNA", "enrichR", "GeneOverlap", "GSEABase", "GSVA"
)

# Function to install packages not yet installed
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# Apply the function to each package
sapply(packages, install_if_missing)

# Bioconductor packages
bioc_packages <- c("SingleCellExperiment", "GSEABase", "GSVA")

# Function to install Bioconductor packages
install_bioc_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(package)
  }
}

# Apply the function to each Bioconductor package
sapply(bioc_packages, install_bioc_if_missing)

# Print message when done
cat("All packages have been checked and installed if necessary.\n")