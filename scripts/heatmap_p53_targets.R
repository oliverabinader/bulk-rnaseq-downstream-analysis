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
output_file <- "results/p53_heatmap.tiff"

# Define p53 target genes
genes <- c(
  "DUSP5","GADD45A","FBXO32","PRDM1","E2F7","CDKN1A","ATF3","PMAIP1",
  "TEP1","CSNK1G1","TP53INP1","CPEB4","BBC3","PLEKHF1","YPEL3",
  "DNAJB2","TRIML2","NUPR1"
)

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

ggsave(
    filename = output_file,
    plot = p,
    device = "tiff",
    width = 6,
    height = 4,
    units = "in",
    dpi = 350,
    compression = "lzw"
  )

# =========================
# Notes:
# - Row scaling highlights relative expression differences
# - Blue = low expression, Red = high expression
# ============================================
