#!/usr/bin/env Rscript
#' Differential Gene Expression Analysis with DESeq2
#' 
#' This script performs comprehensive differential expression analysis using DESeq2
#' on nf-core/rnaseq pipeline output data.
#' 
#' Input: salmon.merged.gene_counts.tsv from nf-core/rnaseq
#' Output: Multiple comparison results, plots, and summary tables
#' 
#' CUSTOMIZATION:
#' To modify sample groups, edit the sample lists below (lines ~35-60) or use:
#' update_sample_groups(new_control = c("sample1", "sample2"), ...)
#' 
#' Author: Generated for CRC differential expression analysis
#' Date: 2025-10-06

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(biomaRt)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(circlize)
  library(readr)
  library(tibble)
  library(stringr)
})

# Set working directory and create output folders
setwd("/Users/kulmanm/KAUST/crc")
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

#' Define sample groups for analysis
#' Modify these lists according to your experimental design
CONTROL_SAMPLES <- c(
  "CRC-10-NT",
  "CRC-13-NT", 
  "CRC-14-NT",
  "CRC-16-NT"
)

TREATMENT_SAMPLES <- c(
  "CRC-10-CT",
  "CRC-13-CT",
  "CRC-14-CT", 
  "CRC-16-CT"
)

# Additional sample groups (organoids)
ORGANOID_PASSAGE_SAMPLES <- c(
  "CRC-10-ENA-P8",
  "CRC-13-ENAS-P3",
  "CRC-14-ENA-P4",
  "CRC-16-ENA-P7"
)

ORGANOID_SPHERE_SAMPLES <- c(
  "CRC-10-ENAS-P8",
  "CRC-13-ENAS-P8",
  "CRC-14-ENAS-P3",
  "CRC-10-ENAS-P8H"
)

#' Helper function to customize sample groups
#' Use this function to easily modify sample assignments before running analysis
#' @param new_control Vector of control sample names
#' @param new_treatment Vector of treatment sample names  
#' @param new_organoid_passage Vector of organoid passage sample names
#' @param new_organoid_sphere Vector of organoid sphere sample names
update_sample_groups <- function(new_control = NULL, new_treatment = NULL, 
                                new_organoid_passage = NULL, new_organoid_sphere = NULL) {
  if (!is.null(new_control)) CONTROL_SAMPLES <<- new_control
  if (!is.null(new_treatment)) TREATMENT_SAMPLES <<- new_treatment
  if (!is.null(new_organoid_passage)) ORGANOID_PASSAGE_SAMPLES <<- new_organoid_passage
  if (!is.null(new_organoid_sphere)) ORGANOID_SPHERE_SAMPLES <<- new_organoid_sphere
  
  cat("Updated sample groups:\n")
  cat("Control:", paste(CONTROL_SAMPLES, collapse = ", "), "\n")
  cat("Treatment:", paste(TREATMENT_SAMPLES, collapse = ", "), "\n")
  cat("Organoid passage:", paste(ORGANOID_PASSAGE_SAMPLES, collapse = ", "), "\n")
  cat("Organoid sphere:", paste(ORGANOID_SPHERE_SAMPLES, collapse = ", "), "\n")
}

#' Print current sample group assignments
print_sample_groups <- function() {
  cat("Current sample group assignments:\n")
  cat("Control samples:", paste(CONTROL_SAMPLES, collapse = ", "), "\n")
  cat("Treatment samples:", paste(TREATMENT_SAMPLES, collapse = ", "), "\n")
  cat("Organoid passage samples:", paste(ORGANOID_PASSAGE_SAMPLES, collapse = ", "), "\n")
  cat("Organoid sphere samples:", paste(ORGANOID_SPHERE_SAMPLES, collapse = ", "), "\n\n")
}

#' Load and prepare count data
#' @param counts_file Path to the gene counts file
#' @param control_samples Vector of control sample names
#' @param treatment_samples Vector of treatment sample names
#' @return List containing counts matrix and sample metadata
load_count_data <- function(counts_file, control_samples = CONTROL_SAMPLES, treatment_samples = TREATMENT_SAMPLES) {
  cat("Loading count data from:", counts_file, "\n")
  
  # Read count data
  counts_raw <- read_tsv(counts_file, show_col_types = FALSE)
  
  # Extract gene info and counts
  gene_info <- counts_raw[, 1:2]  # gene_id and gene_name
  counts_matrix <- counts_raw[, -(1:2)]
  
  # Set row names to gene_id
  rownames(counts_matrix) <- counts_raw$gene_id
  
  # Round counts to integers (DESeq2 requirement)
  counts_matrix <- round(as.matrix(counts_matrix))
  
  # Get all sample names from the count matrix
  all_sample_names <- colnames(counts_matrix)
  
  # Combine all predefined sample groups
  all_predefined_samples <- c(
    control_samples,
    treatment_samples,
    ORGANOID_PASSAGE_SAMPLES,
    ORGANOID_SPHERE_SAMPLES
  )
  
  # Check which predefined samples are actually present in the data
  available_samples <- intersect(all_predefined_samples, all_sample_names)
  missing_samples <- setdiff(all_predefined_samples, all_sample_names)
  
  if (length(missing_samples) > 0) {
    cat("Warning: The following predefined samples are not found in the data:\n")
    cat(paste(missing_samples, collapse = ", "), "\n")
  }
  
  # Use only available samples for analysis
  counts_matrix <- counts_matrix[, available_samples]
  sample_names <- colnames(counts_matrix)
  
  # Create sample metadata based on predefined groups
  sample_metadata <- data.frame(
    sample_id = sample_names,
    stringsAsFactors = FALSE
  )
  
  # Assign conditions based on predefined sample lists
  sample_metadata$condition_simple <- case_when(
    sample_metadata$sample_id %in% control_samples ~ "normal",
    sample_metadata$sample_id %in% treatment_samples ~ "tumor", 
    sample_metadata$sample_id %in% ORGANOID_PASSAGE_SAMPLES ~ "organoid_passage",
    sample_metadata$sample_id %in% ORGANOID_SPHERE_SAMPLES ~ "organoid_sphere",
    TRUE ~ "other"
  )
  
  # Extract patient ID from sample names for batch effect correction
  sample_metadata$patient_id <- str_extract(sample_names, "CRC-\\d+")
  sample_metadata$batch <- sample_metadata$patient_id
  
  # Extract original condition from sample names for reference
  sample_metadata$original_condition <- str_extract(sample_names, "(CT|NT|ENA-P\\d+|ENAS-P\\d+|ENAS-P\\d+H?)$")
  
  rownames(sample_metadata) <- sample_names
  
  cat("Data summary:\n")
  cat("- Genes:", nrow(counts_matrix), "\n")
  cat("- Total samples available:", length(all_sample_names), "\n")
  cat("- Samples used for analysis:", ncol(counts_matrix), "\n")
  cat("- Control samples:", sum(sample_metadata$condition_simple == "normal"), "\n")
  cat("- Treatment samples:", sum(sample_metadata$condition_simple == "tumor"), "\n")
  cat("- Organoid passage samples:", sum(sample_metadata$condition_simple == "organoid_passage"), "\n")
  cat("- Organoid sphere samples:", sum(sample_metadata$condition_simple == "organoid_sphere"), "\n")
  
  # Print sample assignments
  cat("\nSample assignments:\n")
  for (condition in unique(sample_metadata$condition_simple)) {
    samples_in_condition <- sample_metadata$sample_id[sample_metadata$condition_simple == condition]
    cat(paste0(condition, ": ", paste(samples_in_condition, collapse = ", ")), "\n")
  }
  
  return(list(
    counts = counts_matrix,
    metadata = sample_metadata,
    gene_info = gene_info
  ))
}

#' Filter low expression genes
#' @param counts Count matrix
#' @param min_count Minimum count threshold
#' @param min_samples Minimum number of samples
filter_genes <- function(counts, min_count = 10, min_samples = 3) {
  keep <- rowSums(counts >= min_count) >= min_samples
  cat("Filtering genes: keeping", sum(keep), "out of", nrow(counts), "genes\n")
  return(counts[keep, ])
}

#' Create DESeq2 dataset
#' @param counts Filtered count matrix
#' @param metadata Sample metadata
#' @param design Design formula
create_deseq_dataset <- function(counts, metadata, design = ~ batch + condition_simple) {
  cat("Creating DESeq2 dataset with design:", deparse(design), "\n")
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design
  )
  
  # Set reference level (normal as reference for tumor vs normal comparison)
  dds$condition_simple <- relevel(factor(dds$condition_simple), ref = "normal")
  dds$batch <- factor(dds$batch)
  
  return(dds)
}

#' Run DESeq2 analysis
#' @param dds DESeq2 dataset
run_deseq_analysis <- function(dds) {
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds, parallel = TRUE)
  return(dds)
}

#' Extract results for specific contrast
#' @param dds DESeq2 dataset
#' @param contrast Contrast specification
#' @param alpha Significance threshold
extract_deseq_results <- function(dds, contrast, alpha = 0.05, gene_info = NULL) {
  cat("Extracting results for contrast:", paste(contrast, collapse = " vs "), "\n")
  
  res <- results(dds, contrast = contrast, alpha = alpha)
  res_ordered <- res[order(res$padj), ]
  
  # Add gene information if provided
  if (!is.null(gene_info)) {
    res_df <- as.data.frame(res_ordered) %>%
      rownames_to_column("gene_id") %>%
      left_join(gene_info, by = "gene_id")
  } else {
    res_df <- as.data.frame(res_ordered) %>%
      rownames_to_column("gene_id")
  }
  
  # Summary statistics
  cat("Results summary:\n")
  cat("- Total genes:", nrow(res_df), "\n")
  cat("- Significant (padj <", alpha, "):", sum(res_df$padj < alpha, na.rm = TRUE), "\n")
  cat("- Upregulated:", sum(res_df$padj < alpha & res_df$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("- Downregulated:", sum(res_df$padj < alpha & res_df$log2FoldChange < 0, na.rm = TRUE), "\n")
  
  return(res_df)
}

#' Create volcano plot
#' @param results DESeq2 results dataframe
#' @param title Plot title
create_volcano_plot <- function(results, title, alpha = 0.05, fc_threshold = 1) {
  # Prepare data for plotting
  plot_data <- results %>%
    mutate(
      significant = padj < alpha & abs(log2FoldChange) > fc_threshold,
      regulation = case_when(
        padj < alpha & log2FoldChange > fc_threshold ~ "Upregulated",
        padj < alpha & log2FoldChange < -fc_threshold ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Create volcano plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Regulation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Create MA plot
#' @param results DESeq2 results dataframe
#' @param title Plot title
create_ma_plot <- function(results, title, alpha = 0.05) {
  plot_data <- results %>%
    mutate(significant = padj < alpha & !is.na(padj))
  
  p <- ggplot(plot_data, aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "log10(base mean)",
      y = "log2 Fold Change",
      color = paste("Significant (padj <", alpha, ")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Create heatmap of top genes
#' @param dds DESeq2 dataset
#' @param results DESeq2 results
#' @param top_n Number of top genes to include
create_heatmap <- function(dds, results, top_n = 50, title = "Top Differentially Expressed Genes") {
  # Get top genes
  top_genes <- head(results[!is.na(results$padj), ], top_n)$gene_id
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  top_counts <- norm_counts[top_genes, ]
  
  # Log transform and scale
  log_counts <- log2(top_counts + 1)
  scaled_counts <- t(scale(t(log_counts)))
  
  # Create annotation for samples
  annotation_col <- data.frame(
    condition = colData(dds)$condition_simple,
    batch = colData(dds)$batch
  )
  rownames(annotation_col) <- colnames(scaled_counts)
  
  # Create heatmap
  pheatmap(
    scaled_counts,
    annotation_col = annotation_col,
    scale = "none",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title,
    filename = paste0("results/plots/heatmap_", gsub("[^A-Za-z0-9]", "_", title), ".png"),
    width = 12,
    height = 10
  )
}

#' Perform PCA analysis
#' @param dds DESeq2 dataset
perform_pca <- function(dds) {
  # Variance stabilizing transformation
  vst_data <- vst(dds, blind = FALSE)
  
  # PCA plot
  pca_plot <- plotPCA(vst_data, intgroup = c("condition_simple", "batch")) +
    theme_minimal() +
    labs(title = "PCA of Samples") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(pca_plot)
}

#' Main analysis function
main_analysis <- function() {
  cat("=== DESeq2 Differential Expression Analysis ===\n\n")
  
  # Print current sample group assignments
  print_sample_groups()
  
  # Load data
  data <- load_count_data("data/star_salmon/salmon.merged.gene_counts.tsv")
  
  # Filter genes
  filtered_counts <- filter_genes(data$counts)
  
  # Create DESeq2 dataset
  dds <- create_deseq_dataset(filtered_counts, data$metadata)
  
  # Run DESeq2 analysis
  dds <- run_deseq_analysis(dds)
  
  # Perform PCA
  cat("\nCreating PCA plot...\n")
  pca_plot <- perform_pca(dds)
  ggsave("results/plots/pca_analysis.png", pca_plot, width = 10, height = 8, dpi = 300)
  
  # Define comparisons to perform
  comparisons <- list(
    list(contrast = c("condition_simple", "tumor", "normal"), 
         name = "Tumor_vs_Normal"),
    list(contrast = c("condition_simple", "organoid_passage", "normal"), 
         name = "Organoid_Passage_vs_Normal"),
    list(contrast = c("condition_simple", "organoid_sphere", "normal"), 
         name = "Organoid_Sphere_vs_Normal"),
    list(contrast = c("condition_simple", "organoid_sphere", "organoid_passage"), 
         name = "Organoid_Sphere_vs_Passage")
  )
  
  # Perform each comparison
  all_results <- list()
  
  for (comp in comparisons) {
    cat("\n=== Analysis:", comp$name, "===\n")
    
    # Extract results
    results <- extract_deseq_results(dds, comp$contrast, gene_info = data$gene_info)
    all_results[[comp$name]] <- results
    
    # Save results table
    write_csv(results, paste0("results/tables/", comp$name, "_results.csv"))
    
    # Create plots
    volcano_plot <- create_volcano_plot(results, paste("Volcano Plot:", comp$name))
    ma_plot <- create_ma_plot(results, paste("MA Plot:", comp$name))
    
    # Save plots
    ggsave(paste0("results/plots/volcano_", comp$name, ".png"), 
           volcano_plot, width = 10, height = 8, dpi = 300)
    ggsave(paste0("results/plots/ma_", comp$name, ".png"), 
           ma_plot, width = 10, height = 8, dpi = 300)
    
    # Create heatmap for top genes
    create_heatmap(dds, results, title = paste("Top DE Genes:", comp$name))
    
    # Save significant genes
    sig_genes <- results[results$padj < 0.05 & !is.na(results$padj), ]
    write_csv(sig_genes, paste0("results/tables/", comp$name, "_significant_genes.csv"))
  }
  
  # Save DESeq2 object for future use
  saveRDS(dds, "results/deseq2_object.rds")
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results saved in 'results/' directory\n")
  cat("- Tables: results/tables/\n")
  cat("- Plots: results/plots/\n")
  cat("- DESeq2 object: results/deseq2_object.rds\n")
  
  return(list(dds = dds, results = all_results))
}

# Run main analysis if script is executed directly
if (!interactive()) {
  results <- main_analysis()
}
