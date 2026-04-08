# ============================================
# Heatmap of overlapping p53 target genes
# ============================================

# Load required libraries
library(heatmaply)
library(dplyr)

# =========================
# User-defined inputs
# =========================

input_file <- "/path/to/tpm_matrix.csv"
output_file <- "results/p53_heatmap.html"

# Read in the p53 target genes
genes <- readLines("data/gene_sets/p53_targets_reference.txt")

# =========================
# Load data
# =========================

tpm <- read.csv(input_file)

# Subset to target genes
tpm_subset <- tpm %>%
  filter(GeneSymbol %in% genes)

# Set gene names as rownames
rownames(tpm_subset) <- tpm_subset$GeneSymbol

# Remove unnecessary columns if present
tpm_subset <- tpm_subset %>%
  select(-any_of(c("Geneid", "GeneSymbol")))

# Log2 transform
tpm_log <- log2(tpm_subset + 1)

# =========================
# Generate heatmap
# =========================

# OPTION 1: Heatmaply without ggsave
heatmaply(
  as.matrix(tpm_log),
  scale = "row",
  Rowv = TRUE,
  Colv = FALSE,
  showticklabels = c(TRUE, TRUE),
  column_text_angle = 90,
  col = colorRampPalette(c("blue", "white", "red"))(100),
  file = output_file
)

# OTPION 2: Pheatmap with .tiff extension
pheatmap::pheatmap(
  as.matrix(tpm_log),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = "results/p53_heatmap.tiff",
  width = 6,
  height = 8
)

# =========================
# Notes:
# - Row scaling highlights relative expression differences
# - Blue = low expression, Red = high expression
# ============================================
