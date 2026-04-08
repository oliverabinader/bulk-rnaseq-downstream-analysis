# ============================================
# Overlapping Gene Heatmap + Correlation Analysis
# ============================================

# Objective:
# Identify overlapping differentially expressed genes (DEGs) across datasets and visualize:
#   1) Heatmap of overlapping genes
#   2) Correlation of log2 fold changes between conditions

# ============================================
# Load libraries
# ============================================

library(dplyr)
library(readxl)
library(heatmaply)
library(ggplot2)

# ============================================
# User-defined inputs
# ============================================

p53_file <- "data/gene_sets/p53_targets.txt"

de_file <- "data/example_input/DE_results.xlsx"

tpm_file <- "data/example_input/tpm_matrix.csv"

output_dir <- "results/"

# Define sheet names (customize as needed)
sheets <- list(
  CA46_30  = "DEresults.CA46_YX-2-30vsDMSO",
  CA46_228 = "DEresults.CA46_YX-2-228vsDMSO",
  SU6_30   = "DEresults.SU6_YX-2-30vsDMSO",
  SU6_228  = "DEresults.SU6_YX-2-228vsDMSO"
)

# ============================================
# Load data
# ============================================

# Load gene set
p53_targets <- readLines(p53_file)

# Load DE results
de_list <- lapply(sheets, function(sheet) {
  as.data.frame(read_excel(de_file, sheet = sheet))
})

# ============================================
# Filter DEGs
# ============================================

filter_deg <- function(df, fc_thresh = 0.585, fdr_thresh = 0.05) {

  df_up <- df[df$log2FoldChange >= fc_thresh & df$padj < fdr_thresh, ]
  df_dn <- df[df$log2FoldChange <= -fc_thresh & df$padj < fdr_thresh, ]

  return(rbind(df_up, df_dn))
}

deg_filtered <- lapply(de_list, filter_deg)

# ============================================
# Find overlapping genes
# ============================================

# Intersect with gene set
deg_p53 <- lapply(deg_filtered, function(df) {
  intersect(df$GeneSymbol, p53_targets)
})

# Find common genes across all datasets
overlap_genes <- Reduce(intersect, deg_p53)

# ============================================
# Heatmap generation
# ============================================

# Load TPM matrix
tpm <- read.csv(tpm_file)

# Subset to overlapping genes
tpm_subset <- tpm %>%
  filter(GeneSymbol %in% overlap_genes)

rownames(tpm_subset) <- tpm_subset$GeneSymbol

tpm_subset <- tpm_subset %>%
  select(-any_of(c("Geneid", "GeneSymbol")))

# Log transform
tpm_log <- log2(tpm_subset + 1)

# Interactive heatmap
heatmaply(
  as.matrix(tpm_log),
  scale = "row",
  Rowv = TRUE,
  Colv = FALSE,
  showticklabels = c(TRUE, TRUE),
  column_text_angle = 90,
  col = colorRampPalette(c("blue", "white", "red"))(100),
  file = file.path(output_dir, "overlap_heatmap.html")
)

# ============================================
# Correlation analysis
# ============================================

# Example: SU6 comparison between treatments

df1 <- de_list$SU6_30
df2 <- de_list$SU6_228

common_genes_df <- inner_join(
  dplyr::select(df1, Geneid, log2FC_1 = log2FoldChange),
  dplyr::select(df2, Geneid, log2FC_2 = log2FoldChange),
  by = "Geneid"
)

# Spearman correlation
cor_result <- cor.test(
  common_genes_df$log2FC_1,
  common_genes_df$log2FC_2,
  method = "spearman",
  exact = FALSE
)

print(cor_result)

# ============================================
# Plot correlation
# ============================================

p <- ggplot(common_genes_df, aes(x = log2FC_1, y = log2FC_2)) +
  geom_point(color = "orange", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() +
  labs(
    x = "Log2FC Condition 1",
    y = "Log2FC Condition 2",
    title = "Spearman Correlation of Gene Expression"
  )

ggsave(
  filename = file.path(output_dir, "correlation_plot.tiff"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

# ============================================
# Notes:
# - Overlap ensures robust gene selection
# - Correlation assesses consistency between conditions
# ============================================
