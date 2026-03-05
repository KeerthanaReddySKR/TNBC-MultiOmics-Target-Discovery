# ============================================================
# STEP 5 — Differential Expression (limma) + Volcano
# Project: TNBC TCGA-BRCA
# Cohort: expr_1218 (frozen)
# Contrast: TNBC vs non-TNBC
# Data: log2 RSEM (UCSC Xena)
# Purpose:
#   - limma DE with eBayes
#   - Full ranked DE table
#   - Thresholded gene lists
#   - GSEA ranked list
#   - Publication-grade volcano (PDF + PNG)
#   - Summary report
# ============================================================

rm(list = ls())
cat("🧬 STEP 5: limma DE + Volcano — Publication-grade module starting...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(data.table)
  library(matrixStats)
  library(scales)
  library(ggrepel)
})

# -----------------------------
# Paths
# -----------------------------
in_file <- "D:/TNBC_Project_CLEAN/data_processed/expr_1218.rds"

out_dir <- "D:/TNBC_Project_CLEAN/results/DE_TNBC_vs_nonTNBC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Global Figure System (same philosophy as STEP 4)
# -----------------------------

FIG_WIDTH  <- 7.5
FIG_HEIGHT <- 6.0
FIG_DPI    <- 300
BASE_SIZE  <- 14

COLORS <- list(
  DE = c(
    "Up"       = "#D55E00",  # vermillion
    "Down"     = "#0072B2",  # blue
    "NotSig"   = "grey70"
  )
)

theme_project <- function(base_size = 14) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle   = element_text(size = base_size, hjust = 0),
      axis.title      = element_text(face = "bold"),
      axis.text       = element_text(color = "black"),
      legend.title    = element_text(face = "bold"),
      legend.key      = element_blank(),
      panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.6),
      plot.margin     = margin(10, 12, 10, 10)
    )
}

save_figure <- function(plot, filename_base, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI) {
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
    dpi = dpi
  )
}

# -----------------------------
# Load data
# -----------------------------
cat("📥 Loading frozen cohort: expr_1218...\n")
obj <- readRDS(in_file)

expr     <- obj$expr      # genes x samples (log2 RSEM)
labels   <- obj$labels
clinical <- obj$clinical

cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

# -----------------------------
# Sanity checks
# -----------------------------
stopifnot(all(colnames(expr) == labels$sampleID))
stopifnot("is_TNBC" %in% colnames(labels))

table_group <- table(labels$is_TNBC)
cat("Group sizes (FALSE=non-TNBC, TRUE=TNBC):\n")
print(table_group)

if (any(is.na(labels$is_TNBC))) {
  stop("❌ NA values found in is_TNBC labels. This should not happen in frozen cohort.")
}

# -----------------------------
# Remove zero-variance genes (IMPORTANT FIX)
# -----------------------------
cat("🧹 Removing zero-variance genes...\n")

vars <- apply(expr, 1, var)
zero_var_genes <- sum(vars == 0)

cat("Zero-variance genes removed:", zero_var_genes, "\n")

expr <- expr[vars > 0, ]

cat("Genes remaining after filtering:", nrow(expr), "\n")

# -----------------------------
# Design matrix
# -----------------------------
group <- factor(ifelse(labels$is_TNBC, "TNBC", "NonTNBC"),
                levels = c("NonTNBC", "TNBC"))

design <- model.matrix(~ 0 + group)
colnames(design) <- c("NonTNBC", "TNBC")

cat("Design matrix columns:\n")
print(colnames(design))

# -----------------------------
# limma fit
# -----------------------------
cat("🧪 Fitting limma model...\n")

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  TNBC_vs_NonTNBC = TNBC - NonTNBC,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# Extract results
# -----------------------------
cat("📤 Extracting full DE table...\n")

res <- topTable(
  fit2,
  coef = "TNBC_vs_NonTNBC",
  number = Inf,
  sort.by = "P"
)

# Add gene column
res$Gene <- rownames(res)

# Reorder columns
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# Save full table
full_table_file <- file.path(out_dir, "DE_full_TNBC_vs_nonTNBC.tsv")
fwrite(res, full_table_file, sep = "\t")

# -----------------------------
# Thresholds
# -----------------------------
FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1.0

res$DE_class <- "NotSig"
res$DE_class[res$adj.P.Val < FDR_CUTOFF & res$logFC >=  LOGFC_CUTOFF] <- "Up"
res$DE_class[res$adj.P.Val < FDR_CUTOFF & res$logFC <= -LOGFC_CUTOFF] <- "Down"

# -----------------------------
# Save gene lists
# -----------------------------
up_genes   <- res$Gene[res$DE_class == "Up"]
down_genes <- res$Gene[res$DE_class == "Down"]

fwrite(data.table(Gene = up_genes),
       file.path(out_dir, "DE_up_TNBC.tsv"), sep = "\t")

fwrite(data.table(Gene = down_genes),
       file.path(out_dir, "DE_down_TNBC.tsv"), sep = "\t")

# -----------------------------
# Save ranked list for GSEA (FIXED — NO data.table syntax)
# -----------------------------
ranked <- data.frame(
  Gene = res$Gene,
  t = res$t
)

ranked_file <- file.path(out_dir, "DE_ranked_list_for_GSEA.rnk")
fwrite(ranked, ranked_file, sep = "\t", col.names = FALSE)

# -----------------------------
# Summary report
# -----------------------------
qc_report <- c()
qc_report <- c(qc_report, "STEP 5 — Differential Expression (limma)")
qc_report <- c(qc_report, paste0("Date: ", Sys.time()))
qc_report <- c(qc_report, paste0("Genes tested: ", nrow(res)))
qc_report <- c(qc_report, paste0("TNBC samples: ", sum(group == "TNBC")))
qc_report <- c(qc_report, paste0("Non-TNBC samples: ", sum(group == "NonTNBC")))
qc_report <- c(qc_report, paste0("Zero-variance genes removed: ", zero_var_genes))
qc_report <- c(qc_report, paste0("FDR cutoff: ", FDR_CUTOFF))
qc_report <- c(qc_report, paste0("logFC cutoff: ", LOGFC_CUTOFF))
qc_report <- c(qc_report, paste0("Upregulated genes: ", length(up_genes)))
qc_report <- c(qc_report, paste0("Downregulated genes: ", length(down_genes)))

qc_file <- file.path(out_dir, "DE_QC_summary.txt")
writeLines(qc_report, qc_file)

# -----------------------------
# Volcano plot data
# -----------------------------
vol <- res
vol$negLog10FDR <- -log10(vol$adj.P.Val)

# Cap infinite values
vol$negLog10FDR[!is.finite(vol$negLog10FDR)] <- max(vol$negLog10FDR[is.finite(vol$negLog10FDR)], na.rm = TRUE)

# -----------------------------
# Choose top genes to label
# -----------------------------
top_up   <- head(vol[vol$DE_class == "Up",   ][order(vol$adj.P.Val[vol$DE_class == "Up"]), "Gene"], 10)
top_down <- head(vol[vol$DE_class == "Down", ][order(vol$adj.P.Val[vol$DE_class == "Down"]), "Gene"], 10)

label_genes <- unique(c(top_up, top_down))

vol$Label <- ifelse(vol$Gene %in% label_genes, vol$Gene, NA)

# -----------------------------
# Volcano plot
# -----------------------------
p_volcano <- ggplot(vol, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = DE_class), size = 1.6, alpha = 0.85) +
  scale_color_manual(values = COLORS$DE) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3.5,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    na.rm = TRUE
  ) +
  labs(
    title = "Differential expression: TNBC vs non-TNBC",
    subtitle = "limma (eBayes), UCSC Xena log2 RSEM • Threshold: FDR < 0.05 & |log2FC| ≥ 1",
    x = "log2 Fold Change (TNBC / non-TNBC)",
    y = "-log10(FDR)",
    color = "Status"
  ) +
  theme_project(BASE_SIZE)

out_base_volcano <- file.path(out_dir, "Volcano_TNBC_vs_nonTNBC")
save_figure(p_volcano, out_base_volcano)

# -----------------------------
# Final messages
# -----------------------------
cat("\n🏁 STEP 5 COMPLETED SUCCESSFULLY\n")
cat("Results folder:\n", out_dir, "\n")
cat("Full DE table:\n", full_table_file, "\n")
cat("Up genes:", length(up_genes), " Down genes:", length(down_genes), "\n")
cat("Volcano:\n", out_base_volcano, "(.pdf/.png)\n")
cat("GSEA ranked list:\n", ranked_file, "\n")
cat("QC report:\n", qc_file, "\n")
