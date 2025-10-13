#!/usr/bin/env Rscript
#' Manual Package Installation for DESeq2 Analysis
#' 
#' This script provides a step-by-step installation process
#' Run this if the automated setup script fails.

cat("=== Manual DESeq2 Package Installation ===\n\n")

# Set CRAN mirror
cat("Setting CRAN mirror...\n")
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install packages one by one with error handling
install_package_safely <- function(pkg, from = "CRAN") {
  cat("\n--- Installing", pkg, "---\n")
  
  # Check if already installed
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "already installed and loaded successfully\n")
    return(TRUE)
  }
  
  # Try to install
  tryCatch({
    if (from == "CRAN") {
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    } else if (from == "Bioconductor") {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.rstudio.com/")
      }
      BiocManager::install(pkg, dependencies = TRUE)
    }
    
    # Test loading
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("✓", pkg, "installed and loaded successfully\n")
      return(TRUE)
    } else {
      cat("✗", pkg, "installed but failed to load\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("✗ Failed to install", pkg, ":\n")
    cat("  Error:", e$message, "\n")
    cat("  You may need to install this package manually\n")
    return(FALSE)
  })
}

# Core packages (install these first)
cat("\n=== Installing Core Packages ===\n")
core_packages <- c("dplyr", "ggplot2", "readr", "tibble")
core_success <- sapply(core_packages, install_package_safely, "CRAN")

# BiocManager
cat("\n=== Installing BiocManager ===\n")
biocmanager_success <- install_package_safely("BiocManager", "CRAN")

# DESeq2 (most important)
cat("\n=== Installing DESeq2 ===\n")
deseq2_success <- install_package_safely("DESeq2", "Bioconductor")

# Additional CRAN packages
cat("\n=== Installing Additional CRAN Packages ===\n")
additional_cran <- c("stringr", "tidyr", "pheatmap", "RColorBrewer")
cran_success <- sapply(additional_cran, install_package_safely, "CRAN")

# Additional Bioconductor packages
cat("\n=== Installing Additional Bioconductor Packages ===\n")
additional_bioc <- c("clusterProfiler", "org.Hs.eg.db")
bioc_success <- sapply(additional_bioc, install_package_safely, "Bioconductor")

# Optional packages (for advanced features)
cat("\n=== Installing Optional Packages ===\n")
optional_packages <- c("plotly", "DT", "VennDiagram", "UpSetR")
optional_success <- sapply(optional_packages, function(pkg) {
  cat("Attempting to install optional package:", pkg, "\n")
  install_package_safely(pkg, "CRAN")
})

# Summary
cat("\n=== Installation Summary ===\n")
cat("Core packages:", sum(core_success), "/", length(core_success), "successful\n")
cat("BiocManager:", ifelse(biocmanager_success, "✓", "✗"), "\n")
cat("DESeq2:", ifelse(deseq2_success, "✓", "✗"), "\n")
cat("Additional CRAN:", sum(cran_success), "/", length(cran_success), "successful\n")
cat("Additional Bioconductor:", sum(bioc_success), "/", length(bioc_success), "successful\n")
cat("Optional packages:", sum(optional_success), "/", length(optional_success), "successful\n")

# Check minimum requirements
if (deseq2_success && sum(core_success) >= 3) {
  cat("\n✓ Minimum requirements met! You can run the basic DESeq2 analysis.\n")
  
  # Create directory structure
  cat("\nCreating directory structure...\n")
  dir.create("results", showWarnings = FALSE)
  dir.create("results/plots", showWarnings = FALSE)
  dir.create("results/tables", showWarnings = FALSE)
  cat("✓ Directories created\n")
  
  cat("\nYou can now run: Rscript deseq2_analysis.R\n")
} else {
  cat("\n⚠ Minimum requirements not met.\n")
  cat("Please ensure DESeq2 and core packages are installed before proceeding.\n")
  
  if (!deseq2_success) {
    cat("\nTo install DESeq2 manually, run in R:\n")
    cat("if (!require('BiocManager', quietly = TRUE))\n")
    cat("    install.packages('BiocManager')\n")
    cat("BiocManager::install('DESeq2')\n")
  }
}

cat("\n=== Setup Complete ===\n")
