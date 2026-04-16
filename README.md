# Bulk RNA-seq Downstream Analysis

**Author:** Oliver Abinader


## Overview

This repository contains modular downstream analysis workflows for bulk RNA-seq data.  
The focus is on gene-level integration, visualization, and cross-condition comparison.

Analyses include:
- Gene set–based heatmaps (e.g., p53 target genes)
- Robust Rank Aggregation (RRA) across multiple conditions
- Identification of overlapping differentially expressed genes (DEGs)
- Correlation analysis between treatment conditions
- A Venn diagram-based module has also been added to identify overlapping DEGs across conditions

## Analysis Workflow

1. Differential expression results are imported from multiple conditions  
2. Genes are filtered based on:
   - log2 fold-change thresholds  
   - adjusted p-value (FDR) cutoff  
3. Overlapping gene sets (e.g., p53 target genes) are identified  
4. Expression matrices are subset and visualized using heatmaps  
5. Rank-based integration (RRA) is applied across datasets  
6. Correlation analysis is performed between conditions to assess concordance
7. Identification of overlapping DEGs using Venn diagram analysis across conditions

## 📊 Analyses Included

### 1. Gene Set Heatmap Analysis

Visualizes expression patterns of predefined gene sets.

**Features:**
- Uses TPM or normalized expression values
- Log2 transformation applied
- Row-wise scaling for visualization
- Interactive heatmap output (HTML)

### 2. Robust Rank Aggregation (RRA)

Integrates multiple differential expression datasets to identify genes with consistent regulation.

**Identifies:**
- Consistently upregulated genes  
- Consistently downregulated genes  

**Approach:**
- Genes are ranked by log2 fold-change within each dataset
- RRA integrates ranked lists across conditions
- Outputs statistically ranked gene lists

### 3. Overlapping p53 Target Gene Analysis

Identifies significantly differentially expressed p53-associated genes shared across conditions.

**Steps:**
- Filter DEGs using:
  - |log2FC| ≥ 0.585 (≈ 1.5-fold change)
  - FDR < 0.05
- Intersect DEGs with curated p53 gene set
- Identify genes shared across multiple datasets
- Generate heatmaps of overlapping genes

### 4. Correlation Analysis Between Conditions

Assesses transcriptional concordance between treatment conditions.

**Features:**
- Spearman correlation of log2 fold-changes
- Gene-level comparison between conditions
- Scatter plot visualization with regression line

### 5. Overlapping DEG analysis (Venn diagram)

**Steps:**
- Differential expression results are used directly
- Genes are filtered using:
  - log2 fold-change thresholds
  - adjusted p-value (FDR cutoff)
- Upregulated and downregulated genes are separated
- Overlaps across conditions are identified using Venn diagrams
- Shared gene sets are exported for downstream analysis

## 📥 Input Data Requirements

### Differential Expression (DE) Files
Each dataset must contain:

- `GeneSymbol`
- `Geneid` (recommended)
- `log2FoldChange`
- `padj` (adjusted p-value)

### Expression Matrix (TPM or normalized counts)

- Rows = genes
- Columns = samples/conditions
- Must include `GeneSymbol` column

### Gene Sets

Example:
- p53 target gene list (`.txt` file)


## 📊 Outputs
-  Heatmaps (HTML / TIFF)
-  RRA ranked gene tables (CSV)
-  Correlation plots (TIFF/PNG)
-  Overlapping gene lists (TXT)
-  Venn diagrams (PNG) and shared gene lists (Excel)

## Installation

install.packages(c("heatmaply", "RobustRankAggreg", "dplyr"))


## 📂 Repository Structure

bulk-rnaseq-downstream-analysis/
- scripts/
   - RRA_analysis.R
   - heatmap_p53_targets.R
   - p53_overlap_and_correlation.R
   - Venn_DEG_analysis.R 
- gene_sets/
   - p53_targets_reference.txt
   - p53_targets_overlap_deg.txt
- results/
   - heatmaps/
   - correlation/
   - rra/
- README.md
