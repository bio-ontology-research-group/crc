#!/usr/bin/env Rscript
#' Setup Script for DESeq2 Analysis Pipeline
#' 
#' This script installs all required R packages and sets up the environment
#' for running the DESeq2 differential expression analysis.
#' 
#' Run this script before running the main analysis:
#' Rscript setup_deseq2_environment.R

cat("=== DESeq2 Analysis Environment Setup ===\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))
cat("CRAN mirror set to:", getOption("repos")["CRAN"], "\n\n")

# Function to install packages if not already installed
install_if_missing <- function(packages, from = "CRAN") {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "from", from, "...\n")
      tryCatch({
        if (from == "CRAN") {
          install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
        } else if (from == "Bioconductor") {
          BiocManager::install(pkg, dependencies = TRUE)
        }
        cat("✓", pkg, "installed successfully\n")
      }, error = function(e) {
        cat("✗ Failed to install", pkg, ":", e$message, "\n")
      })
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Install BiocManager if not available
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  tryCatch({
    install.packages("BiocManager", repos = "https://cran.rstudio.com/")
    cat("✓ BiocManager installed successfully\n")
  }, error = function(e) {
    cat("✗ Failed to install BiocManager:", e$message, "\n")
    cat("Please install BiocManager manually: install.packages('BiocManager')\n")
  })
} else {
  cat("✓ BiocManager already available\n")
}

# CRAN packages
cat("\n--- Installing CRAN packages ---\n")
cran_packages <- c(
  "dplyr",
  "ggplot2", 
  "readr",
  "tibble",
  "stringr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "VennDiagram",
  "UpSetR",
  "plotly",
  "DT",
  "htmlwidgets",
  "ggrepel"
)

install_if_missing(cran_packages, "CRAN")

# Bioconductor packages
cat("\n--- Installing Bioconductor packages ---\n")
bioc_packages <- c(
  "DESeq2",
  "biomaRt",
  "EnhancedVolcano",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot"
)

install_if_missing(bioc_packages, "Bioconductor")

# Load all packages to test installation
cat("\n--- Testing package loading ---\n")
all_packages <- c(cran_packages, bioc_packages)
failed_packages <- c()

for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    failed_packages <- c(failed_packages, pkg)
  }
}

# Report results
cat("\n=== Installation Summary ===\n")
if (length(failed_packages) == 0) {
  cat("✓ All packages installed successfully!\n")
  cat("You can now run the DESeq2 analysis scripts.\n")
} else {
  cat("⚠ The following packages failed to install:\n")
  for (pkg in failed_packages) {
    cat("  -", pkg, "\n")
  }
  cat("\nPlease install these packages manually before proceeding.\n")
}

# Create results directory structure
cat("\n--- Setting up directory structure ---\n")
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/gene_sets", showWarnings = FALSE)

cat("✓ Directory structure created\n")

# Check data file
cat("\n--- Checking data files ---\n")
if (file.exists("data/star_salmon/salmon.merged.gene_counts.tsv")) {
  cat("✓ Gene counts file found\n")
} else {
  cat("⚠ Gene counts file not found at: data/star_salmon/salmon.merged.gene_counts.tsv\n")
  cat("Please ensure your data file is in the correct location.\n")
}

cat("\n=== Setup Complete ===\n")
cat("Next steps:\n")
cat("1. Run main analysis: Rscript deseq2_analysis.R\n")
cat("2. Run advanced analysis: Rscript deseq2_advanced_analysis.R\n")
cat("3. Run utilities: Rscript deseq2_utilities.R\n")
cat("\nSee DESeq2_Analysis_README.md for detailed instructions.\n")
