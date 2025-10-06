#!/usr/bin/env Rscript
#' Advanced DESeq2 Analysis and Functional Enrichment
#' 
#' This script provides additional analysis functions including:
#' - Gene set enrichment analysis
#' - Pathway analysis
#' - Time series analysis for organoid passages
#' - Custom visualization functions
#' 
#' Author: Generated for CRC differential expression analysis
#' Date: 2025-10-06

# Load additional libraries for advanced analysis
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggrepel)
  library(VennDiagram)
  library(UpSetR)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Vector of gene symbols
#' @return Vector of Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  gene_map <- bitr(gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)
  return(gene_map)
}

#' Perform GO enrichment analysis
#' @param gene_list Vector of gene symbols
#' @param universe Background gene universe
#' @param ont GO ontology (BP, MF, CC)
perform_go_enrichment <- function(gene_list, universe = NULL, ont = "BP", title = "") {
  cat("Performing GO enrichment analysis for", length(gene_list), "genes\n")
  
  # Convert to Entrez IDs
  gene_entrez <- convert_to_entrez(gene_list)
  
  if (nrow(gene_entrez) == 0) {
    cat("Warning: No genes could be converted to Entrez IDs\n")
    return(NULL)
  }
  
  # Convert universe if provided
  if (!is.null(universe)) {
    universe_entrez <- convert_to_entrez(universe)$ENTREZID
  } else {
    universe_entrez <- NULL
  }
  
  # Perform GO enrichment
  ego <- enrichGO(gene = gene_entrez$ENTREZID,
                  universe = universe_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    cat("No significant GO terms found\n")
    return(NULL)
  }
  
  cat("Found", nrow(ego@result), "significant GO terms\n")
  return(ego)
}

#' Perform KEGG pathway enrichment
#' @param gene_list Vector of gene symbols
#' @param universe Background gene universe
perform_kegg_enrichment <- function(gene_list, universe = NULL, title = "") {
  cat("Performing KEGG pathway enrichment for", length(gene_list), "genes\n")
  
  # Convert to Entrez IDs
  gene_entrez <- convert_to_entrez(gene_list)
  
  if (nrow(gene_entrez) == 0) {
    cat("Warning: No genes could be converted to Entrez IDs\n")
    return(NULL)
  }
  
  # Convert universe if provided
  if (!is.null(universe)) {
    universe_entrez <- convert_to_entrez(universe)$ENTREZID
  } else {
    universe_entrez <- NULL
  }
  
  # Perform KEGG enrichment
  kegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                     universe = universe_entrez,
                     organism = 'hsa',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  
  if (is.null(kegg) || nrow(kegg@result) == 0) {
    cat("No significant KEGG pathways found\n")
    return(NULL)
  }
  
  cat("Found", nrow(kegg@result), "significant KEGG pathways\n")
  return(kegg)
}

#' Create enrichment plots
#' @param enrichment_result clusterProfiler result object
#' @param title Plot title
#' @param output_prefix Output file prefix
create_enrichment_plots <- function(enrichment_result, title, output_prefix) {
  if (is.null(enrichment_result)) return(NULL)
  
  # Dot plot
  p1 <- dotplot(enrichment_result, showCategory = 20) + 
    ggtitle(paste("Dot Plot:", title)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Bar plot
  p2 <- barplot(enrichment_result, showCategory = 20) + 
    ggtitle(paste("Bar Plot:", title)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Gene-concept network (if terms available)
  if (nrow(enrichment_result@result) >= 5) {
    p3 <- cnetplot(enrichment_result, categorySize = "pvalue", foldChange = NULL) +
      ggtitle(paste("Network Plot:", title)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0("results/plots/", output_prefix, "_network.png"), 
           p3, width = 12, height = 10, dpi = 300)
  }
  
  # Save plots
  ggsave(paste0("results/plots/", output_prefix, "_dotplot.png"), 
         p1, width = 10, height = 8, dpi = 300)
  ggsave(paste0("results/plots/", output_prefix, "_barplot.png"), 
         p2, width = 10, height = 8, dpi = 300)
  
  return(list(dotplot = p1, barplot = p2))
}

#' Analyze overlap between gene sets
#' @param gene_lists Named list of gene vectors
#' @param title Plot title
analyze_gene_overlap <- function(gene_lists, title = "Gene Set Overlap") {
  cat("Analyzing overlap between", length(gene_lists), "gene sets\n")
  
  # Create UpSet plot
  upset_data <- fromList(gene_lists)
  
  png(paste0("results/plots/upset_", gsub("[^A-Za-z0-9]", "_", title), ".png"), 
      width = 12, height = 8, units = "in", res = 300)
  upset(upset_data, 
        sets = names(gene_lists),
        order.by = "freq",
        keep.order = TRUE,
        main.bar.color = "steelblue",
        sets.bar.color = "darkred")
  dev.off()
  
  # Create Venn diagram if 2-3 sets
  if (length(gene_lists) %in% 2:3) {
    venn_colors <- c("red", "blue", "green")[1:length(gene_lists)]
    
    venn.diagram(
      x = gene_lists,
      category.names = names(gene_lists),
      filename = paste0("results/plots/venn_", gsub("[^A-Za-z0-9]", "_", title), ".png"),
      output = TRUE,
      imagetype = "png",
      height = 2000,
      width = 2000,
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      col = venn_colors,
      fill = alpha(venn_colors, 0.3),
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans"
    )
  }
  
  # Calculate pairwise overlaps
  overlap_matrix <- matrix(0, nrow = length(gene_lists), ncol = length(gene_lists),
                          dimnames = list(names(gene_lists), names(gene_lists)))
  
  for (i in 1:length(gene_lists)) {
    for (j in 1:length(gene_lists)) {
      if (i != j) {
        overlap <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
        overlap_matrix[i, j] <- overlap
      } else {
        overlap_matrix[i, j] <- length(gene_lists[[i]])
      }
    }
  }
  
  return(overlap_matrix)
}

#' Create expression heatmap across conditions
#' @param dds DESeq2 object
#' @param genes Vector of genes to plot
#' @param title Plot title
create_expression_heatmap <- function(dds, genes, title = "Gene Expression Heatmap") {
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Subset to genes of interest
  plot_data <- norm_counts[intersect(genes, rownames(norm_counts)), ]
  
  if (nrow(plot_data) == 0) {
    cat("Warning: No genes found in dataset\n")
    return(NULL)
  }
  
  # Log transform and scale
  log_data <- log2(plot_data + 1)
  scaled_data <- t(scale(t(log_data)))
  
  # Create sample annotation
  sample_annotation <- data.frame(
    condition = colData(dds)$condition_simple,
    patient = colData(dds)$batch,
    row.names = colnames(scaled_data)
  )
  
  # Create heatmap
  pheatmap(
    scaled_data,
    annotation_col = sample_annotation,
    scale = "none",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = title,
    filename = paste0("results/plots/expression_heatmap_", gsub("[^A-Za-z0-9]", "_", title), ".png"),
    width = 12,
    height = max(6, nrow(scaled_data) * 0.3)
  )
}

#' Analyze organoid progression (passage analysis)
#' @param dds DESeq2 object
#' @param results_list List of DESeq2 results
analyze_organoid_progression <- function(dds, results_list) {
  cat("Analyzing organoid progression patterns...\n")
  
  # Get genes that change across passages
  passage_genes <- results_list[["Organoid_Passage_vs_Normal"]]
  sphere_genes <- results_list[["Organoid_Sphere_vs_Normal"]]
  sphere_vs_passage <- results_list[["Organoid_Sphere_vs_Passage"]]
  
  # Find genes with significant changes
  sig_passage <- passage_genes[passage_genes$padj < 0.05 & !is.na(passage_genes$padj), ]
  sig_sphere <- sphere_genes[sphere_genes$padj < 0.05 & !is.na(sphere_genes$padj), ]
  sig_progression <- sphere_vs_passage[sphere_vs_passage$padj < 0.05 & !is.na(sphere_vs_passage$padj), ]
  
  # Identify progression patterns
  # Genes that increase with progression
  increasing_genes <- sig_progression[sig_progression$log2FoldChange > 0, ]$gene_id
  # Genes that decrease with progression  
  decreasing_genes <- sig_progression[sig_progression$log2FoldChange < 0, ]$gene_id
  
  # Create expression heatmap for progression genes
  if (length(increasing_genes) > 0) {
    create_expression_heatmap(dds, head(increasing_genes, 50), "Increasing_with_Organoid_Progression")
  }
  
  if (length(decreasing_genes) > 0) {
    create_expression_heatmap(dds, head(decreasing_genes, 50), "Decreasing_with_Organoid_Progression")
  }
  
  # Save progression gene lists
  write_csv(data.frame(gene_id = increasing_genes), "results/tables/organoid_increasing_genes.csv")
  write_csv(data.frame(gene_id = decreasing_genes), "results/tables/organoid_decreasing_genes.csv")
  
  return(list(
    increasing = increasing_genes,
    decreasing = decreasing_genes
  ))
}

#' Main advanced analysis function
run_advanced_analysis <- function() {
  cat("=== Advanced DESeq2 Analysis ===\n\n")
  
  # Load results from main analysis
  if (!file.exists("results/deseq2_object.rds")) {
    cat("Error: Please run the main DESeq2 analysis first\n")
    return(NULL)
  }
  
  dds <- readRDS("results/deseq2_object.rds")
  
  # Load results tables
  result_files <- list.files("results/tables", pattern = "_results.csv", full.names = TRUE)
  results_list <- list()
  
  for (file in result_files) {
    name <- str_remove(basename(file), "_results.csv")
    results_list[[name]] <- read_csv(file, show_col_types = FALSE)
  }
  
  # Get background gene universe
  all_genes <- rownames(dds)
  
  # Perform enrichment analysis for each comparison
  for (comp_name in names(results_list)) {
    cat("\n=== Enrichment Analysis:", comp_name, "===\n")
    
    results <- results_list[[comp_name]]
    sig_genes <- results[results$padj < 0.05 & !is.na(results$padj), ]
    
    if (nrow(sig_genes) > 0) {
      # Get upregulated and downregulated genes
      up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]$gene_name
      down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]$gene_name
      
      # Remove NA values
      up_genes <- up_genes[!is.na(up_genes)]
      down_genes <- down_genes[!is.na(down_genes)]
      
      if (length(up_genes) > 5) {
        # GO enrichment for upregulated genes
        go_up <- perform_go_enrichment(up_genes, all_genes, "BP", paste(comp_name, "Upregulated"))
        if (!is.null(go_up)) {
          create_enrichment_plots(go_up, paste(comp_name, "UP GO-BP"), paste0(comp_name, "_up_GO"))
          write_csv(go_up@result, paste0("results/tables/", comp_name, "_upregulated_GO.csv"))
        }
        
        # KEGG enrichment for upregulated genes
        kegg_up <- perform_kegg_enrichment(up_genes, all_genes, paste(comp_name, "Upregulated"))
        if (!is.null(kegg_up)) {
          create_enrichment_plots(kegg_up, paste(comp_name, "UP KEGG"), paste0(comp_name, "_up_KEGG"))
          write_csv(kegg_up@result, paste0("results/tables/", comp_name, "_upregulated_KEGG.csv"))
        }
      }
      
      if (length(down_genes) > 5) {
        # GO enrichment for downregulated genes
        go_down <- perform_go_enrichment(down_genes, all_genes, "BP", paste(comp_name, "Downregulated"))
        if (!is.null(go_down)) {
          create_enrichment_plots(go_down, paste(comp_name, "DOWN GO-BP"), paste0(comp_name, "_down_GO"))
          write_csv(go_down@result, paste0("results/tables/", comp_name, "_downregulated_GO.csv"))
        }
        
        # KEGG enrichment for downregulated genes
        kegg_down <- perform_kegg_enrichment(down_genes, all_genes, "BP", paste(comp_name, "Downregulated"))
        if (!is.null(kegg_down)) {
          create_enrichment_plots(kegg_down, paste(comp_name, "DOWN KEGG"), paste0(comp_name, "_down_KEGG"))
          write_csv(kegg_down@result, paste0("results/tables/", comp_name, "_downregulated_KEGG.csv"))
        }
      }
    }
  }
  
  # Analyze gene overlaps between comparisons
  cat("\n=== Gene Overlap Analysis ===\n")
  sig_gene_lists <- list()
  for (comp_name in names(results_list)) {
    results <- results_list[[comp_name]]
    sig_genes <- results[results$padj < 0.05 & !is.na(results$padj), ]$gene_name
    sig_genes <- sig_genes[!is.na(sig_genes)]
    if (length(sig_genes) > 0) {
      sig_gene_lists[[comp_name]] <- sig_genes
    }
  }
  
  if (length(sig_gene_lists) > 1) {
    overlap_matrix <- analyze_gene_overlap(sig_gene_lists, "Significant_Genes_All_Comparisons")
    write.csv(overlap_matrix, "results/tables/gene_overlap_matrix.csv")
  }
  
  # Organoid progression analysis
  if ("Organoid_Sphere_vs_Passage" %in% names(results_list)) {
    cat("\n=== Organoid Progression Analysis ===\n")
    progression_genes <- analyze_organoid_progression(dds, results_list)
    
    # Enrichment analysis for progression genes
    if (length(progression_genes$increasing) > 5) {
      go_inc <- perform_go_enrichment(progression_genes$increasing, all_genes, "BP", "Organoid Increasing")
      if (!is.null(go_inc)) {
        create_enrichment_plots(go_inc, "Organoid Progression Increasing GO-BP", "organoid_increasing_GO")
        write_csv(go_inc@result, "results/tables/organoid_increasing_GO.csv")
      }
    }
    
    if (length(progression_genes$decreasing) > 5) {
      go_dec <- perform_go_enrichment(progression_genes$decreasing, all_genes, "BP", "Organoid Decreasing")
      if (!is.null(go_dec)) {
        create_enrichment_plots(go_dec, "Organoid Progression Decreasing GO-BP", "organoid_decreasing_GO")
        write_csv(go_dec@result, "results/tables/organoid_decreasing_GO.csv")
      }
    }
  }
  
  cat("\n=== Advanced Analysis Complete ===\n")
  cat("Additional results saved in 'results/' directory\n")
}

# Run advanced analysis if script is executed directly
if (!interactive()) {
  run_advanced_analysis()
}
