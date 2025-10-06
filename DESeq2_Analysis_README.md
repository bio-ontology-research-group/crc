# DESeq2 Differential Expression Analysis for CRC Data

This repository contains comprehensive R scripts for differential gene expression analysis using DESeq2 on nf-core/rnaseq pipeline output data.

## Overview

The analysis pipeline consists of four main R scripts:

1. **`deseq2_analysis.R`** - Main differential expression analysis
2. **`deseq2_advanced_analysis.R`** - Functional enrichment and pathway analysis
3. **`deseq2_utilities.R`** - Utility functions and interactive visualizations
4. **Custom analysis examples** - Templates for specific comparisons

## Data Structure

Your data should follow this structure:
```
data/star_salmon/salmon.merged.gene_counts.tsv
```

The gene counts file should contain:
- Column 1: `gene_id` - Ensembl gene IDs
- Column 2: `gene_name` - Gene symbols
- Columns 3+: Sample count data

## Sample Naming Convention

Samples are automatically parsed based on the naming pattern: `CRC-##-CONDITION`

Where:
- `CRC-##` = Patient ID (e.g., CRC-10, CRC-13, CRC-14)
- `CONDITION` = Sample condition:
  - `CT` = Tumor tissue
  - `NT` = Normal tissue  
  - `ENA-P#` = Organoid passage
  - `ENAS-P#` = Organoid sphere
  - `ENAS-P#H` = Organoid sphere variant

## Installation

### Required R Packages

Install the required packages before running the analysis:

```r
# CRAN packages
install.packages(c(
  "dplyr", "ggplot2", "readr", "tibble", "stringr",
  "pheatmap", "RColorBrewer", "VennDiagram", "UpSetR",
  "plotly", "DT", "htmlwidgets"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "biomaRt", "EnhancedVolcano", "ComplexHeatmap",
  "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot"
))
```

## Usage

### 1. Main Differential Expression Analysis

Run the primary DESeq2 analysis:

```bash
cd /Users/kulmanm/KAUST/crc
Rscript deseq2_analysis.R
```

This script will:
- Load and filter gene count data
- Create DESeq2 dataset with batch correction
- Perform differential expression analysis for multiple comparisons:
  - Tumor vs Normal
  - Organoid Passage vs Normal
  - Organoid Sphere vs Normal
  - Organoid Sphere vs Passage
- Generate volcano plots, MA plots, and heatmaps
- Create PCA analysis
- Save results tables and plots

**Output:**
- `results/tables/` - CSV files with differential expression results
- `results/plots/` - Visualization plots (PNG format)
- `results/deseq2_object.rds` - Saved DESeq2 object for further analysis

### 2. Advanced Functional Analysis

Run enrichment and pathway analysis:

```bash
Rscript deseq2_advanced_analysis.R
```

This script performs:
- GO (Gene Ontology) enrichment analysis
- KEGG pathway enrichment
- Gene set overlap analysis
- Organoid progression analysis
- Creates network plots and enrichment visualizations

**Output:**
- `results/tables/*_GO.csv` - GO enrichment results
- `results/tables/*_KEGG.csv` - KEGG pathway results
- `results/plots/*_network.png` - Gene-pathway network plots
- `results/plots/upset_*.png` - Gene set overlap plots

### 3. Utility Analysis and Interactive Visualizations

Generate additional utilities and interactive plots:

```bash
Rscript deseq2_utilities.R
```

This script creates:
- Quality control plots
- Interactive volcano plots (HTML)
- Interactive results tables (HTML)
- Gene set exports for external tools
- Summary reports

**Output:**
- `results/plots/qc_*.png` - Quality control plots
- `results/plots/interactive_*.html` - Interactive visualizations
- `results/gene_sets/*.txt` - Gene lists for external analysis
- `results/analysis_summary.csv` - Overall analysis summary

## Results Structure

After running all scripts, your results directory will contain:

```
results/
├── tables/
│   ├── Tumor_vs_Normal_results.csv
│   ├── Organoid_Passage_vs_Normal_results.csv
│   ├── Organoid_Sphere_vs_Normal_results.csv
│   ├── Organoid_Sphere_vs_Passage_results.csv
│   ├── *_significant_genes.csv
│   ├── *_GO.csv
│   ├── *_KEGG.csv
│   ├── interactive_*.html
│   └── analysis_summary.csv
├── plots/
│   ├── volcano_*.png
│   ├── ma_*.png
│   ├── heatmap_*.png
│   ├── pca_analysis.png
│   ├── qc_*.png
│   ├── *_network.png
│   ├── upset_*.png
│   └── interactive_*.html
├── gene_sets/
│   ├── *_all_significant.txt
│   ├── *_upregulated.txt
│   ├── *_downregulated.txt
│   └── *_top100.txt
└── deseq2_object.rds
```

## Key Features

### Statistical Analysis
- **Batch correction**: Accounts for patient-to-patient variation
- **Multiple comparisons**: Automatically performs relevant contrasts
- **FDR correction**: Uses Benjamini-Hochberg adjustment
- **Robust filtering**: Removes low-expression genes

### Visualizations
- **Volcano plots**: Show fold change vs significance
- **MA plots**: Display mean vs fold change relationship
- **Heatmaps**: Hierarchical clustering of top genes
- **PCA plots**: Sample clustering and batch effects
- **Interactive plots**: HTML-based exploration tools

### Functional Analysis
- **GO enrichment**: Biological process, molecular function, cellular component
- **KEGG pathways**: Metabolic and signaling pathways
- **Gene overlap**: Venn diagrams and UpSet plots
- **Network visualization**: Gene-pathway relationships

## Customization

### Custom Contrasts

To perform specific comparisons, use the utility functions:

```r
# Load analysis data
source("deseq2_utilities.R")
data <- load_analysis_data()
dds <- data$dds

# Custom contrast example: Compare specific patients
custom_results <- perform_custom_contrast(
  dds, 
  contrast = c("batch", "CRC-14", "CRC-10"),
  contrast_name = "Patient_14_vs_10"
)
```

### Modify Analysis Parameters

Edit the main script to adjust:

```r
# In deseq2_analysis.R
# Change significance thresholds
alpha <- 0.01  # More stringent p-value
fc_threshold <- 1.5  # Higher fold-change requirement

# Modify design formula
design <- ~ condition_simple  # Remove batch correction
design <- ~ batch + treatment + condition_simple  # Add treatment factor
```

### Add New Comparisons

Add new contrasts to the comparisons list:

```r
# In deseq2_analysis.R
comparisons <- list(
  # Existing comparisons...
  list(contrast = c("condition_simple", "tumor", "organoid_sphere"), 
       name = "Tumor_vs_Organoid_Sphere")
)
```

## Interpretation Guide

### Results Tables
- **baseMean**: Average normalized count across all samples
- **log2FoldChange**: Log2 fold change (positive = upregulated in first condition)
- **lfcSE**: Standard error of log2 fold change
- **stat**: Test statistic
- **pvalue**: Raw p-value
- **padj**: Adjusted p-value (FDR)

### Significance Thresholds
- **Significant genes**: padj < 0.05
- **Highly significant**: padj < 0.01
- **Biologically relevant**: |log2FoldChange| > 1 (2-fold change)

## Troubleshooting

### Common Issues

1. **Missing packages**: Install all required packages from CRAN and Bioconductor
2. **Memory issues**: Reduce number of genes or samples if needed
3. **Convergence errors**: May indicate too few replicates or extreme outliers

### Performance Tips

- **Parallel processing**: DESeq2 uses multiple cores automatically
- **Memory management**: Large datasets may require more RAM
- **File paths**: Use absolute paths for reliability

## Support

For questions about the analysis pipeline:
1. Check the R script comments for parameter explanations
2. Refer to the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
3. Review the [clusterProfiler documentation](https://yulab-smu.top/biomedical-knowledge-mining-book/) for enrichment analysis

## Citation

If you use these scripts in your research, please cite:

- Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550.
- Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284-287.

## License

This analysis pipeline is provided as-is for research purposes. Please ensure proper attribution when using or modifying these scripts.
