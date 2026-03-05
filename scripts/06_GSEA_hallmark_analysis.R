# ============================================================
# STEP 6: GSEA (Hallmark) on TNBC vs non-TNBC ranked list
# Publication-grade module
# ============================================================

rm(list = ls())
cat("🧬 Starting STEP 6: Hallmark GSEA...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# -----------------------------
# Paths
# -----------------------------
ranked_file <- "D:/TNBC_Project_CLEAN/results/DE_TNBC_vs_nonTNBC/DE_ranked_list_for_GSEA.rnk"
out_dir     <- "D:/TNBC_Project_CLEAN/results/GSEA_Hallmark"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load ranked list
# -----------------------------
cat("📥 Loading ranked list...\n")
rnk <- fread(ranked_file, header = FALSE)
colnames(rnk) <- c("Gene", "Score")

# Remove duplicated genes (keep highest absolute score)
rnk <- rnk[order(-abs(Score))]
rnk <- rnk[!duplicated(Gene)]

# Build named vector
ranks <- rnk$Score
names(ranks) <- rnk$Gene

cat("Genes in ranked list:", length(ranks), "\n")

# -----------------------------
# Load MSigDB Hallmark gene sets
# -----------------------------
cat("📚 Loading MSigDB Hallmark gene sets...\n")

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_sets <- split(msig_h$gene_symbol, msig_h$gs_name)

cat("Number of Hallmark pathways:", length(hallmark_sets), "\n")

# -----------------------------
# Run fgsea
# -----------------------------
cat("🚀 Running fgsea...\n")

set.seed(123)

fg <- fgsea(
  pathways = hallmark_sets,
  stats    = ranks,
  nperm    = 10000,
  minSize  = 15,
  maxSize  = 500
)

fg <- fg[order(fg$padj), ]

# Convert to data.table
fg_dt <- as.data.table(fg)

# -----------------------------
# Save full results
# -----------------------------
full_file <- file.path(out_dir, "GSEA_Hallmark_full_results.tsv")
fwrite(fg_dt, full_file, sep = "\t")

cat("Saved full GSEA results to:", full_file, "\n")

# -----------------------------
# Split upregulated vs downregulated
# -----------------------------
fg_up   <- fg_dt[NES > 0 & padj < 0.05]
fg_down <- fg_dt[NES < 0 & padj < 0.05]

# Sort by NES
fg_up   <- fg_up[order(-NES)]
fg_down <- fg_down[order(NES)]

# Save top tables
fwrite(fg_up,   file.path(out_dir, "GSEA_Hallmark_UP.tsv"),   sep = "\t")
fwrite(fg_down, file.path(out_dir, "GSEA_Hallmark_DOWN.tsv"), sep = "\t")

# -----------------------------
# QC Summary
# -----------------------------
qc_file <- file.path(out_dir, "GSEA_QC_summary.txt")

qc_lines <- c(
  "STEP 6: Hallmark GSEA QC Summary",
  "===============================",
  paste("Date:", Sys.time()),
  "",
  paste("Genes in ranked list:", length(ranks)),
  paste("Hallmark pathways tested:", length(hallmark_sets)),
  "",
  paste("Significant pathways (FDR < 0.05):", sum(fg_dt$padj < 0.05)),
  paste("  Upregulated in TNBC (NES > 0):", nrow(fg_up)),
  paste("  Downregulated in TNBC (NES < 0):", nrow(fg_down)),
  "",
  "Top 5 UP pathways:",
  paste0("  - ", fg_up$pathway[1:min(5, nrow(fg_up))]),
  "",
  "Top 5 DOWN pathways:",
  paste0("  - ", fg_down$pathway[1:min(5, nrow(fg_down))])
)

writeLines(qc_lines, qc_file)

cat("Saved QC summary to:", qc_file, "\n")

# -----------------------------
# Plotting theme (paper style)
# -----------------------------
theme_paper <- function() {
  theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# -----------------------------
# Function to save PDF + PNG
# -----------------------------
save_both <- function(plot, filename_base, width = 7, height = 6) {
  ggsave(
    filename = paste0(filename_base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    filename = paste0(filename_base, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# -----------------------------
# Barplot of top pathways
# -----------------------------
cat("🎨 Generating enrichment barplots...\n")

# Take top 10 up and down
topN <- 10

plot_df <- rbind(
  fg_up[1:min(topN, nrow(fg_up))],
  fg_down[1:min(topN, nrow(fg_down))]
)

plot_df$pathway <- str_replace(plot_df$pathway, "HALLMARK_", "")
plot_df$pathway <- str_replace_all(plot_df$pathway, "_", " ")

plot_df$Direction <- ifelse(plot_df$NES > 0, "Up in TNBC", "Down in TNBC")

p <- ggplot(plot_df, aes(x = reorder(pathway, NES), y = NES, fill = Direction)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Up in TNBC" = "#D55E00", "Down in TNBC" = "#0072B2")) +
  labs(
    title = "Hallmark GSEA: TNBC vs non-TNBC",
    subtitle = "Top enriched pathways (FDR < 0.05)",
    x = "",
    y = "Normalized Enrichment Score (NES)",
    fill = ""
  ) +
  theme_paper()

out_base <- file.path(out_dir, "GSEA_Hallmark_TopPathways")
save_both(p, out_base)

cat("Saved GSEA barplot to:", out_base, "(.pdf/.png)\n")

# -----------------------------
# Finish
# -----------------------------
cat("✅ STEP 6 Hallmark GSEA completed successfully.\n")
