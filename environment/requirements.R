# R Package Requirements for Multiome 10x Analysis
# Run this script to install all required packages

# ============================================================================
# CRAN Packages
# ============================================================================

cran_packages <- c(
  # Data manipulation
  "tidyverse",
  "dplyr",
  "purrr",
  "tibble",
  "readr",
  
  # Visualization
  "ggplot2",
  "patchwork",
  "cowplot",
  "RColorBrewer",
  "viridis",
  "pheatmap",
  "ComplexHeatmap",
  
  # Statistical analysis
  "broom",
  "nlme",
  
  # Utilities
  "devtools",
  "BiocManager"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    cat("✓ Installed:", pkg, "\n")
  } else {
    cat("✓ Already installed:", pkg, "\n")
  }
}

# ============================================================================
# Bioconductor Packages
# ============================================================================

bioc_packages <- c(
  # Single-cell analysis
  "SingleCellExperiment",
  "scater",
  "scran",
  
  # ATAC-seq
  "GenomicRanges",
  "IRanges",
  "chromVAR",
  "motifmatchr",
  "TFBSTools",
  "JASPAR2020",
  
  # Annotations
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "org.Hs.eg.db",
  
  # Utilities
  "BiocGenerics",
  "S4Vectors"
)

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
    cat("✓ Installed:", pkg, "\n")
  } else {
    cat("✓ Already installed:", pkg, "\n")
  }
}

# ============================================================================
# GitHub Packages
# ============================================================================

cat("\nInstalling packages from GitHub...\n")

# Seurat (if not already installed from CRAN)
if (!require("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

# Signac for ATAC-seq analysis
if (!require("Signac", quietly = TRUE)) {
  devtools::install_github("stuart-lab/signac")
}

# SeuratWrappers for integration methods
if (!require("SeuratWrappers", quietly = TRUE)) {
  devtools::install_github("satijalab/seurat-wrappers")
}

# ============================================================================
# Verify Installation
# ============================================================================

cat("\n" , rep("=", 60), "\n", sep = "")
cat("VERIFYING INSTALLATION\n")
cat(rep("=", 60), "\n", sep = "")

key_packages <- c("Seurat", "Signac", "tidyverse", "chromVAR", "GenomicRanges")

all_installed <- TRUE
for (pkg in key_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "- FAILED\n")
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\n✓ All key packages installed successfully!\n")
} else {
  cat("\n✗ Some packages failed to install. Please check errors above.\n")
}

# ============================================================================
# Session Info
# ============================================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SESSION INFO\n")
cat(rep("=", 60), "\n", sep = "")
print(sessionInfo())