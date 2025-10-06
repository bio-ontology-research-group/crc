#!/usr/bin/env Rscript
#' DESeq2 Utility Functions and Custom Analysis
#' 
#' This script provides utility functions for:
#' - Custom contrasts and comparisons
#' - Time series analysis
#' - Quality control plots
#' - Data export functions
#' - Interactive analysis helpers
#' 
#' Author: Generated for CRC differential expression analysis
#' Date: 2025-10-06

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(stringr)
  library(plotly)
  library(DT)
  library(htmlwidgets)
})

#' Load DESeq2 object and results
#' @return List containing DESeq2 object and results
load_analysis_data <- function() {
  if (!file.exists("results/deseq2_object.rds")) {
    stop("Please run the main DESeq2 analysis first")
  }
  
  dds <- readRDS("results/deseq2_object.rds")
  
  # Load all results
  result_files <- list.files("results/tables", pattern = "_results.csv", full.names = TRUE)
  results_list <- list()
  
  for (file in result_files) {
    name <- str_remove(basename(file), "_results.csv")
    results_list[[name]] <- read_csv(file, show_col_types = FALSE)
  }
  
  return(list(dds = dds, results = results_list))
}

#' Perform custom contrast analysis
#' @param dds DESeq2 object
#' @param contrast_vector Custom contrast vector
#' @param contrast_name Name for the contrast
#' @param alpha Significance threshold
perform_custom_contrast <- function(dds, contrast_vector, contrast_name, alpha = 0.05) {
  cat("Performing custom contrast:", contrast_name, "\n")
  
  # Extract results
  res <- results(dds, contrast = contrast_vector, alpha = alpha)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)
  
  # Add gene names if available
  if ("gene_name" %in% colnames(mcols(dds))) {
    gene_info <- data.frame(
      gene_id = rownames(dds),
      gene_name = mcols(dds)$gene_name
    )
    res_df <- res_df %>% left_join(gene_info, by = "gene_id")
  }
  
  # Save results
  write_csv(res_df, paste0("results/tables/custom_", contrast_name, "_results.csv"))
  
  # Summary
  cat("Results summary for", contrast_name, ":\n")
  cat("- Total genes:", nrow(res_df), "\n")
  cat("- Significant (padj <", alpha, "):", sum(res_df$padj < alpha, na.rm = TRUE), "\n")
  cat("- Upregulated:", sum(res_df$padj < alpha & res_df$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("- Downregulated:", sum(res_df$padj < alpha & res_df$log2FoldChange < 0, na.rm = TRUE), "\n")
  
  return(res_df)
}

#' Create interactive volcano plot
#' @param results DESeq2 results dataframe
#' @param title Plot title
#' @param alpha Significance threshold
#' @param fc_threshold Fold change threshold
create_interactive_volcano <- function(results, title, alpha = 0.05, fc_threshold = 1) {
  # Prepare data
  plot_data <- results %>%
    mutate(
      significant = padj < alpha & abs(log2FoldChange) > fc_threshold,
      regulation = case_when(
        padj < alpha & log2FoldChange > fc_threshold ~ "Upregulated",
        padj < alpha & log2FoldChange < -fc_threshold ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      # Create hover text
      hover_text = paste0(
        "Gene: ", ifelse(!is.na(gene_name), gene_name, gene_id), "<br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "padj: ", formatC(padj, format = "e", digits = 2), "<br>",
        "baseMean: ", round(baseMean, 2)
      )
    ) %>%
    filter(!is.na(padj) & !is.na(log2FoldChange))
  
  # Create plotly volcano plot
  p <- plot_ly(
    plot_data,
    x = ~log2FoldChange,
    y = ~-log10(padj),
    color = ~regulation,
    colors = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray"),
    text = ~hover_text,
    hovertemplate = "%{text}<extra></extra>",
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.7)
  ) %>%
    layout(
      title = title,
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10(adjusted p-value)"),
      hovermode = "closest"
    ) %>%
    add_lines(
      x = c(-fc_threshold, -fc_threshold), 
      y = c(0, max(-log10(plot_data$padj), na.rm = TRUE)),
      line = list(dash = "dash", color = "black", width = 1),
      showlegend = FALSE, inherit = FALSE
    ) %>%
    add_lines(
      x = c(fc_threshold, fc_threshold), 
      y = c(0, max(-log10(plot_data$padj), na.rm = TRUE)),
      line = list(dash = "dash", color = "black", width = 1),
      showlegend = FALSE, inherit = FALSE
    ) %>%
    add_lines(
      x = c(min(plot_data$log2FoldChange, na.rm = TRUE), max(plot_data$log2FoldChange, na.rm = TRUE)),
      y = c(-log10(alpha), -log10(alpha)),
      line = list(dash = "dash", color = "black", width = 1),
      showlegend = FALSE, inherit = FALSE
    )
  
  # Save interactive plot
  output_file <- paste0("results/plots/interactive_volcano_", gsub("[^A-Za-z0-9]", "_", title), ".html")
  saveWidget(p, output_file, selfcontained = TRUE)
  
  return(p)
}

#' Create interactive results table
#' @param results DESeq2 results dataframe
#' @param title Table title
create_interactive_table <- function(results, title) {
  # Prepare data for display
  display_data <- results %>%
    mutate(
      log2FoldChange = round(log2FoldChange, 3),
      lfcSE = round(lfcSE, 3),
      stat = round(stat, 3),
      pvalue = formatC(pvalue, format = "e", digits = 2),
      padj = formatC(padj, format = "e", digits = 2),
      baseMean = round(baseMean, 2)
    ) %>%
    select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    filter(!is.na(padj))
  
  # Create DT table
  dt_table <- datatable(
    display_data,
    caption = title,
    filter = "top",
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      autoWidth = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),
    extensions = 'Buttons'
  ) %>%
    formatStyle(
      "padj",
      backgroundColor = styleInterval(c(0.01, 0.05), c("lightgreen", "lightyellow", "white"))
    ) %>%
    formatStyle(
      "log2FoldChange",
      backgroundColor = styleInterval(c(-1, 1), c("lightblue", "white", "lightcoral"))
    )
  
  # Save interactive table
  output_file <- paste0("results/tables/interactive_", gsub("[^A-Za-z0-9]", "_", title), ".html")
  saveWidget(dt_table, output_file, selfcontained = TRUE)
  
  return(dt_table)
}

#' Generate quality control plots
#' @param dds DESeq2 object
generate_qc_plots <- function(dds) {
  cat("Generating quality control plots...\n")
  
  # Sample distances
  vst_data <- vst(dds, blind = FALSE)
  sample_dists <- dist(t(assay(vst_data)))
  sample_dist_matrix <- as.matrix(sample_dists)
  
  # Distance heatmap
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    main = "Sample Distance Heatmap",
    filename = "results/plots/qc_sample_distances.png",
    width = 10,
    height = 8
  )
  
  # Dispersion plot
  png("results/plots/qc_dispersion.png", width = 10, height = 8, units = "in", res = 300)
  plotDispEsts(dds, main = "Dispersion Estimates")
  dev.off()
  
  # MA plot for all genes
  png("results/plots/qc_ma_plot.png", width = 10, height = 8, units = "in", res = 300)
  plotMA(dds, main = "MA Plot - All Genes")
  dev.off()
  
  # Count distribution
  count_data <- as.data.frame(counts(dds, normalized = TRUE)) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
    mutate(log_count = log10(count + 1))
  
  p_dist <- ggplot(count_data, aes(x = log_count, color = sample)) +
    geom_density(alpha = 0.7) +
    labs(
      title = "Distribution of Normalized Counts",
      x = "log10(normalized count + 1)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("results/plots/qc_count_distribution.png", p_dist, width = 10, height = 6, dpi = 300)
  
  cat("QC plots saved in results/plots/\n")
}

#' Export gene sets for external analysis
#' @param results_list List of DESeq2 results
#' @param alpha Significance threshold
export_gene_sets <- function(results_list, alpha = 0.05) {
  cat("Exporting gene sets for external analysis...\n")
  
  # Create directory for gene sets
  dir.create("results/gene_sets", showWarnings = FALSE)
  
  for (comp_name in names(results_list)) {
    results <- results_list[[comp_name]]
    sig_genes <- results[results$padj < alpha & !is.na(results$padj), ]
    
    if (nrow(sig_genes) > 0) {
      # All significant genes
      write_lines(sig_genes$gene_name[!is.na(sig_genes$gene_name)], 
                  paste0("results/gene_sets/", comp_name, "_all_significant.txt"))
      
      # Upregulated genes
      up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
      write_lines(up_genes$gene_name[!is.na(up_genes$gene_name)], 
                  paste0("results/gene_sets/", comp_name, "_upregulated.txt"))
      
      # Downregulated genes
      down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
      write_lines(down_genes$gene_name[!is.na(down_genes$gene_name)], 
                  paste0("results/gene_sets/", comp_name, "_downregulated.txt"))
      
      # Top 100 genes by significance
      top_genes <- head(sig_genes, 100)
      write_lines(top_genes$gene_name[!is.na(top_genes$gene_name)], 
                  paste0("results/gene_sets/", comp_name, "_top100.txt"))
    }
  }
  
  cat("Gene sets exported to results/gene_sets/\n")
}

#' Generate summary report
#' @param results_list List of DESeq2 results
#' @param alpha Significance threshold
generate_summary_report <- function(results_list, alpha = 0.05) {
  cat("Generating summary report...\n")
  
  # Create summary table
  summary_data <- data.frame()
  
  for (comp_name in names(results_list)) {
    results <- results_list[[comp_name]]
    
    total_genes <- nrow(results)
    sig_genes <- sum(results$padj < alpha, na.rm = TRUE)
    up_genes <- sum(results$padj < alpha & results$log2FoldChange > 0, na.rm = TRUE)
    down_genes <- sum(results$padj < alpha & results$log2FoldChange < 0, na.rm = TRUE)
    
    summary_data <- rbind(summary_data, data.frame(
      Comparison = comp_name,
      Total_Genes = total_genes,
      Significant_Genes = sig_genes,
      Upregulated = up_genes,
      Downregulated = down_genes,
      Percent_Significant = round(sig_genes / total_genes * 100, 2)
    ))
  }
  
  # Save summary
  write_csv(summary_data, "results/analysis_summary.csv")
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  print(summary_data)
  
  return(summary_data)
}

#' Main utility analysis function
run_utility_analysis <- function() {
  cat("=== DESeq2 Utility Analysis ===\n\n")
  
  # Load data
  data <- load_analysis_data()
  dds <- data$dds
  results_list <- data$results
  
  # Generate QC plots
  generate_qc_plots(dds)
  
  # Create interactive plots for main comparisons
  for (comp_name in names(results_list)) {
    cat("Creating interactive visualizations for:", comp_name, "\n")
    
    # Interactive volcano plot
    create_interactive_volcano(results_list[[comp_name]], comp_name)
    
    # Interactive results table
    create_interactive_table(results_list[[comp_name]], comp_name)
  }
  
  # Export gene sets
  export_gene_sets(results_list)
  
  # Generate summary report
  summary_data <- generate_summary_report(results_list)
  
  cat("\n=== Utility Analysis Complete ===\n")
  cat("Interactive visualizations and utilities saved in results/\n")
  
  return(list(summary = summary_data, dds = dds, results = results_list))
}

#' Custom analysis examples
#' 
#' Example usage for custom contrasts:
#' 
#' # Load data
#' data <- load_analysis_data()
#' dds <- data$dds
#' 
#' # Perform custom contrast (example: specific patient comparison)
#' custom_contrast <- c("batch", "CRC-14", "CRC-10")  # Compare patient 14 vs 10
#' results <- perform_custom_contrast(dds, custom_contrast, "Patient_14_vs_10")
#' 
#' # Create specific subset analysis
#' # Filter samples to specific conditions
#' tumor_samples <- dds$condition_simple == "tumor"
#' dds_tumor <- dds[, tumor_samples]
#' dds_tumor$batch <- droplevels(dds_tumor$batch)
#' dds_tumor <- DESeq(dds_tumor)
#' 
#' # Compare tumors from different patients
#' tumor_results <- perform_custom_contrast(dds_tumor, c("batch", "CRC-14", "CRC-10"), "Tumor_Patient14_vs_10")

# Run utility analysis if script is executed directly
if (!interactive()) {
  run_utility_analysis()
}
