# Bulk RNA-seq Downstream Analysis

**Author:** Oliver Abinader  

## Overview
This repository contains downstream analyses for bulk RNA-seq data, focusing on:
- Heatmap visualization of selected gene sets (e.g., p53 target genes)
- Robust Rank Aggregation (RRA) to identify genes consistently upregulated or downregulated across multiple conditions

These analyses are designed to be modular and adaptable to different datasets.


## 📊 Analyses Included

### 1. Heatmap of Target Genes
- Visualizes expression of predefined gene sets
- Uses TPM (or normalized counts)
- Applies log transformation and row scaling
- Generates an interactive heatmap

### 2. Robust Rank Aggregation (RRA)
- Integrates multiple differential expression results
- Identifies genes consistently:
  - Upregulated
  - Downregulated
- Based on ranked log2 fold changes

### 3. Overlapping p53 Target Gene Heatmap
- Identifies significantly differentially expressed genes (DEGs)
- Filters based on:
  - log2FC threshold (e.g., ±0.585 ≈ 1.5-fold)
  - adjusted p-value (FDR < 0.05)
- Intersects DEGs across datasets and predefined gene sets
- Generates heatmaps for overlapping genes

### 4. Correlation Analysis Between Conditions
- Computes Spearman correlation between log2 fold changes
- Identifies concordance between treatments
- Visualizes results using scatter plots with regression line

## Input Requirements

### Heatmap Input
- CSV/TSV file with:
  - `GeneSymbol` column
  - Expression values (TPM or normalized counts)

### RRA Input
- Differential expression files containing:
  - `GeneSymbol`
  - `log2FoldChange`


## Installation

install.packages(c("heatmaply", "RobustRankAggreg", "dplyr"))


## Usage
- Run Heatmap: source("scripts/heatmap_p53_targets.R")
- Run RRA Analysis: source("scripts/RRA_analysis.R")


## ⚠️ Notes
- Update file paths in scripts before running
- Input data is not included due to size and/or confidentiality
