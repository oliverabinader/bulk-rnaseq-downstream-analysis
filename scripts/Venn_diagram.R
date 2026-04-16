# ============================================================
# Differential Expression Venn Diagram Pipeline
# Author: Oliver Abinader
# ============================================================

# ---------------------------
# Load libraries
# ---------------------------
library(dplyr)
library(ggVennDiagram)
library(readr)
library(writexl)

# ---------------------------
# USER SETTINGS
# ---------------------------

output_dir <- "/path/to/data_visualization"

fc_cutoff <- 0.585
padj_cutoff <- 0.05

gene_col <- "GeneID"

# ---------------------------
# INPUT FILES
# ---------------------------

# This section defines the differential expression (DE) result files used for Venn diagram analysis.
# The naming scheme encodes the experimental design:
# NP vs AF → tissue/type or group comparison
# DBM vs BW → treatment/condition comparison

files <- list(
  NP_DBM = "data/DEG.NP_DBMvsControl.tsv",
  NP_BW  = "data/DEG.NP_BWvsControl.tsv",
  AF_DBM = "data/DEG.AF_DBMvsControl.tsv",
  AF_BW  = "data/DEG.AF_BWvsControl.tsv"
)

# ---------------------------
# LOAD DATA
# ---------------------------
data_list <- lapply(files, function(f) {
  read.csv(f, check.names = FALSE)
})

# ---------------------------
# FUNCTION: extract DE genes
# ---------------------------
get_genes <- function(df, direction = c("up", "down")) {

  direction <- match.arg(direction)

  if (direction == "up") {
    df %>%
      filter(log2FoldChange >= fc_cutoff & padj < padj_cutoff) %>%
      pull(all_of(gene_col))
  } else {
    df %>%
      filter(log2FoldChange <= -fc_cutoff & padj < padj_cutoff) %>%
      pull(all_of(gene_col))
  }
}

# ---------------------------
# EXTRACT UP & DOWN GENES
# ---------------------------

up_list <- lapply(data_list, get_genes, direction = "up")
down_list <- lapply(data_list, get_genes, direction = "down")

# Clean names
names(up_list) <- names(files)
names(down_list) <- names(files)

# ---------------------------
# FUNCTION: make venn plot
# ---------------------------
make_venn <- function(gene_list, title) {

  ggVennDiagram(
    gene_list,
    label = "count",
    label_alpha = 0,
    label_size = 7
  ) +
    scale_fill_gradient(low = "white", high = "coral") +
    ggtitle(title) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank()
    )
}

# ---------------------------
# CREATE VENN PLOTS
# ---------------------------

venn_up <- make_venn(up_list, "Upregulated Gene Overlap (All Conditions)")
venn_dn <- make_venn(down_list, "Downregulated Gene Overlap (All Conditions)")

# ---------------------------
# SAVE PLOTS
# ---------------------------

ggsave(
  filename = file.path(output_dir, "Venn_UP.png"),
  plot = venn_up,
  device = "tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 350,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "Venn_DOWN.png"),
  plot = venn_dn,
  device = "tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 350,
  compression = "lzw"
)

# ---------------------------
# EXPORT COMMON GENES
# ---------------------------

common_up <- Reduce(intersect, up_list)
common_down <- Reduce(intersect, down_list)

write_xlsx(
  list(
    Upregulated_Common = data.frame(Genes = common_up),
    Downregulated_Common = data.frame(Genes = common_down)
  ),
  file.path(output_dir, "Common_DEG_sets.xlsx")
)

# ---------------------------
# SUMMARY
# ---------------------------

cat("\n==================== SUMMARY ====================\n")

cat("\nUpregulated genes per condition:\n")
print(sapply(up_list, length))

cat("\nDownregulated genes per condition:\n")
print(sapply(down_list, length))

cat("\nCommon upregulated genes:", length(common_up), "\n")
cat("Common downregulated genes:", length(common_down), "\n")

cat("\nOutput saved in:", output_dir, "\n")
