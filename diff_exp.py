#!/usr/bin/env python3
"""
Differential Expression Analysis with PyDESeq2
Author: Your Name
Date: 2025

This script performs differential expression analysis using PyDESeq2,
which is a Python implementation of the popular R DESeq2 package.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('default')
sns.set_palette("husl")

def load_and_prepare_data(counts_file, sample_metadata=None):
    """
    Load count data and prepare sample metadata
    """
    # Load count data from nf-core/rnaseq
    print("Loading count data...")
    treatment_columns = ['CRC-10-ENAS-P8', 'CRC-13-ENAS-P8', 'CRC-14-ENAS-P3']
    control_columns = ['CRC-10-ENA-P8','CRC-13-ENAS-P3', 'CRC-14-ENA-P4']
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    # Subset count data into control and treatment groups
    control_counts = counts_df[control_columns]
    treatment_counts = counts_df[treatment_columns]
    counts_df = pd.concat([control_counts, treatment_counts], axis=1)

    counts_df = counts_df.astype(int)

    # Clean column names (remove everything after first dot)
    counts_df.columns = [col.split('.')[0] for col in counts_df.columns]
    
    # Create sample metadata if not provided
    if sample_metadata is None:
        sample_names = counts_df.columns.tolist()
        
        # assuming first half are controls, second half are treatments
        n_samples = len(sample_names)
        conditions = ['control'] * (n_samples//2) + ['treatment'] * (n_samples - n_samples//2)
        
        sample_metadata = pd.DataFrame({
            'sample': sample_names,
            'condition': conditions,
            'batch': ['batch1'] * n_samples  # Add batch info if needed
        })
        sample_metadata.set_index('sample', inplace=True)
    
    print(f"Loaded {counts_df.shape[0]} genes across {counts_df.shape[1]} samples")
    print("Sample metadata:")
    print(sample_metadata)
    
    return counts_df, sample_metadata

def run_deseq2_analysis(counts_df, sample_metadata, design_formula="~ condition"):
    """
    Run differential expression analysis using PyDESeq2
    """
    print("\nRunning DESeq2 analysis...")
    
    # Filter low count genes (keep genes with at least 10 counts total)
    keep_genes = counts_df.sum(axis=1) >= 10
    counts_filtered = counts_df[keep_genes].T
    print(counts_filtered)
    print(f"Kept {counts_filtered.shape[1]} genes after filtering")
    
    # Create DESeq2 dataset
    inference = DefaultInference(n_cpus=4)  # Adjust CPUs as needed
    
    dds = DeseqDataSet(
        counts=counts_filtered,
        metadata=sample_metadata,
        design_factors=design_formula.replace("~ ", "").split(" + "),
        refit_cooks=True,
        inference=inference
    )
    
    # Run DESeq2
    dds.deseq2()
    
    return dds

def extract_results(dds, contrast=("condition", "treatment", "control"), alpha=0.05):
    """
    Extract differential expression results
    """
    print(f"\nExtracting results for contrast: {contrast}")
    
    stat_res = DeseqStats(dds, contrast=contrast, inference=dds.inference)
    stat_res.summary()
    
    # Get results as DataFrame
    results_df = stat_res.results_df.copy()
    results_df['gene_id'] = results_df.index
    # Filter significant genes
    sig_genes = results_df[
        (results_df['padj'] < alpha) & 
        (abs(results_df['log2FoldChange']) > 1)
    ].copy()
    sig_genes = sig_genes.sort_values('padj')
    
    print(f"Total genes tested: {len(results_df)}")
    print(f"Significantly upregulated genes: {sum(sig_genes['log2FoldChange'] > 1)}")
    print(f"Significantly downregulated genes: {sum(sig_genes['log2FoldChange'] < -1)}")
    
    return results_df, sig_genes


def main():
    """
    Main analysis pipeline
    """
    # File paths - MODIFY THESE PATHS
    counts_file = "data/star_salmon/salmon.merged.gene_counts.tsv"  # Path to your count file
    
    # Optional: load custom sample metadata
    # sample_metadata = pd.read_csv("sample_metadata.csv", index_col=0)
    sample_metadata = None  # Will be auto-generated
    
    try:
        # 1. Load data
        counts_df, sample_metadata = load_and_prepare_data(counts_file, sample_metadata)
        
        # 2. Run DESeq2 analysis
        dds = run_deseq2_analysis(counts_df, sample_metadata)
        
        # 3. Extract results
        results_df, sig_genes = extract_results(dds, 
                                               contrast=("condition", "treatment", "control"))
        
        # 4. Save results
        print("\nSaving results...")
        results_df.to_csv('differential_expression_results.csv')
        sig_genes.to_csv('significant_genes.csv')
        
        print("\n=== Analysis Complete! ===")
        print("Files generated:")
        print("- differential_expression_results.csv: All results")
        print("- significant_genes.csv: Significant genes only")
        
        # Display summary
        print(f"\nSummary:")
        print(f"Total genes: {len(results_df)}")
        print(f"Significant genes: {len(sig_genes)}")
        print(f"Upregulated: {sum(sig_genes['log2FoldChange'] > 1)}")
        print(f"Downregulated: {sum(sig_genes['log2FoldChange'] < -1)}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        print("Please check your file paths and data format.")

if __name__ == "__main__":
    # Install required packages if needed
    print("Installing required packages...")
    import subprocess
    import sys
    
    packages = ['pydeseq2', 'pandas', 'numpy', 'matplotlib', 'seaborn', 'scikit-learn', 'scipy']
    
    for package in packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    
    main()