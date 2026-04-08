# ============================================
# Objective:
# Identify genes that are consistently upregulated or downregulated across multiple RNA-seq comparisons using Robust Rank Aggregation (RRA).
#
# Approach:
# - Genes are ranked within each dataset based on log2 fold-change
# - Separate rankings are generated for:
#     1) Upregulated genes (positive log2FC)
#     2) Downregulated genes (negative log2FC)
# - RRA is applied to identify genes with concordant rankings across all datasets
#
# How RRA works:
# Each dataset is converted into a ranked gene list based on log2FC.
# These ranked lists are stored separately for upregulated and downregulated genes.
#
# The RRA algorithm then evaluates all lists together to identify genes that consistently appear at the top (or bottom) across datasets.
# This does NOT merge datasets directly, but instead integrates their rankings to detect robust, concordant signals.
#
# Output:
# - Ranked gene lists with associated significance scores (FDR)
# ============================================

# Load package
if (!require("RobustRankAggreg")) install.packages("RobustRankAggreg")
library(RobustRankAggreg)

# =========================
# User-defined inputs
# =========================

input_files <- list(
  dataset1 = "/path/to/DE_1.csv",
  dataset2 = "/path/to/DE_2.csv",
  dataset3 = "/path/to/DE_3.csv",
  dataset4 = "/path/to/DE_4.csv"
)
# Each file should contain two columns: GeneSymbol and log2FC

output_dir <- "results/"

# =========================
# Load data
# =========================

data_list <- lapply(input_files, read.csv)

# =========================
# Generate ranked lists
# =========================

rank_lists_pos <- list()
rank_lists_neg <- list()

for (name in names(data_list)) {

  df <- data_list[[name]]

  # -------------------------
  # Positive log2FC genes
  # -------------------------
  pos <- df[df$log2FoldChange > 0, ]
  rank_lists_pos[[name]] <- pos[order(-pos$log2FoldChange), "GeneSymbol"]

  # -------------------------
  # Negative log2FC genes
  # -------------------------
  neg <- df[df$log2FoldChange < 0, ]
  rank_lists_neg[[name]] <- neg[order(neg$log2FoldChange), "GeneSymbol"]
}

# =========================
# Run RRA
# =========================

# Upregulated genes
rra_pos <- aggregateRanks(glist = rank_lists_pos, method = "RRA", full = TRUE)
rra_pos$FDR <- p.adjust(rra_pos$Score, method = "BH", n = nrow(rra_pos)))

# Downregulated genes
rra_neg <- aggregateRanks(glist = rank_lists_neg, method = "RRA", full = TRUE)
rra_neg$FDR <- p.adjust(rra_neg$Score, method = "BH", n = nrow(rra_neg)))

# =========================
# Save results
# =========================

write.csv(rra_pos, file.path(output_dir, "RRA_positive.csv"), row.names = FALSE)
write.csv(rra_neg, file.path(output_dir, "RRA_negative.csv"), row.names = FALSE)

# =========================
# Interpretation
# =========================
# rra_pos:
#   Genes consistently upregulated across datasets
#
# rra_neg:
#   Genes consistently downregulated across datasets
#
# Lower Score / FDR:
#   Higher confidence in consistent ranking
# ============================================
