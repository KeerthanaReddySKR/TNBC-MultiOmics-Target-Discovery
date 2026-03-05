# ============================================================
# STEP 4 — QC & PCA (Publication-grade figure module)
# Project: TNBC TCGA-BRCA
# Cohort: expr_1218
# Purpose:
#   - Global QC
#   - PCA on top variable genes
#   - Publication-quality figures (PDF + PNG)
#   - Variance table
#   - QC summary report
# ============================================================

rm(list = ls())
cat("🧬 STEP 4: QC & PCA — Publication-grade module starting...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(matrixStats)
  library(scales)
  library(grid)
})

# -----------------------------
# Paths
# -----------------------------
in_file <- "D:/TNBC_Project_CLEAN/data_processed/expr_1218.rds"

out_dir <- "D:/TNBC_Project_CLEAN/results/PCA_expr1218"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Global Figure System
# -----------------------------

# ---- Dimensions (Nature-style) ----
FIG_WIDTH  <- 7.5   # inches
FIG_HEIGHT <- 6.0   # inches
FIG_DPI    <- 300

# ---- Base font size ----
BASE_SIZE <- 14

# ---- Color system (Okabe–Ito / Nature-safe) ----
COLORS <- list(
  TNBC = c(
    "TNBC"     = "#D55E00",  # vermillion
    "Non-TNBC" = "#0072B2"   # blue
  ),
  PAM50 = c(
    "Basal"  = "#D55E00",
    "Her2"   = "#CC79A7",
    "LumA"   = "#009E73",
    "LumB"   = "#56B4E9",
    "Normal" = "#999999"
  )
)

# ---- Project-wide ggplot theme ----
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

# ---- Central export function ----
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

expr     <- obj$expr      # genes x samples
labels   <- obj$labels
clinical <- obj$clinical

cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

# -----------------------------
# Sanity checks
# -----------------------------
stopifnot(all(colnames(expr) == labels$sampleID))
stopifnot(all(colnames(expr) == clinical$sampleID))

# -----------------------------
# Basic QC
# -----------------------------
qc_report <- c()
qc_report <- c(qc_report, paste0("STEP 4 QC & PCA REPORT"))
qc_report <- c(qc_report, paste0("Date: ", Sys.time()))
qc_report <- c(qc_report, paste0("Genes: ", nrow(expr)))
qc_report <- c(qc_report, paste0("Samples: ", ncol(expr)))

# Check missing values
na_count <- sum(is.na(expr))
qc_report <- c(qc_report, paste0("Total NA values in expression matrix: ", na_count))

# -----------------------------
# Select top variable genes
# -----------------------------
cat("🔍 Selecting top variable genes for PCA...\n")

gene_var <- rowVars(as.matrix(expr))
names(gene_var) <- rownames(expr)

TOP_N <- 5000

if (length(gene_var) < TOP_N) {
  stop("Not enough genes for PCA.")
}

top_genes <- names(sort(gene_var, decreasing = TRUE))[1:TOP_N]
expr_top <- expr[top_genes, ]

qc_report <- c(qc_report, paste0("Top variable genes used for PCA: ", nrow(expr_top)))

# -----------------------------
# PCA
# -----------------------------
cat("📊 Running PCA (center = TRUE, scale = TRUE)...\n")

pca <- prcomp(t(expr_top), center = TRUE, scale. = TRUE)

# -----------------------------
# Variance explained
# -----------------------------
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

var_df <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance = var_explained,
  Percent = var_explained * 100,
  Cumulative = cumsum(var_explained) * 100
)

var_file <- file.path(out_dir, "PCA_variance_expr1218.tsv")
fwrite(var_df, var_file, sep = "\t")

qc_report <- c(qc_report, "Variance explained (first 5 PCs):")
for (i in 1:5) {
  qc_report <- c(qc_report, paste0(
    "  PC", i, ": ",
    round(var_df$Percent[i], 2), "% (cumulative ",
    round(var_df$Cumulative[i], 2), "%)"
  ))
}

# -----------------------------
# Build PCA data.frame
# -----------------------------
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

# Add TNBC
pca_df$is_TNBC <- labels$is_TNBC[match(pca_df$Sample, labels$sampleID)]
pca_df$TNBC_label <- ifelse(pca_df$is_TNBC, "TNBC", "Non-TNBC")

# Add PAM50
pca_df$PAM50 <- clinical$PAM50Call_RNAseq[match(pca_df$Sample, clinical$sampleID)]

# -----------------------------
# Outlier detection (simple distance-based)
# -----------------------------
d2 <- pca_df$PC1^2 + pca_df$PC2^2
threshold <- quantile(d2, 0.99)
outliers <- pca_df$Sample[d2 > threshold]

qc_report <- c(qc_report, paste0("Outliers detected (99% radius): ", length(outliers)))

outlier_file <- file.path(out_dir, "PCA_outliers_expr1218.txt")
writeLines(outliers, outlier_file)

# -----------------------------
# Axis labels
# -----------------------------
pc1_lab <- paste0("PC1 (", round(var_df$Percent[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_df$Percent[2], 1), "%)")

# -----------------------------
# Plot 1: TNBC
# -----------------------------
p_tnbc <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TNBC_label)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_color_manual(values = COLORS$TNBC) +
  labs(
    title = "PCA of TCGA-BRCA transcriptomes",
    subtitle = "Top 5,000 variable genes • Colored by TNBC status",
    x = pc1_lab,
    y = pc2_lab,
    color = "Subtype"
  ) +
  theme_project(BASE_SIZE)

out_base_tnbc <- file.path(out_dir, "PCA_expr1218_TNBC")
save_figure(p_tnbc, out_base_tnbc)

# -----------------------------
# Plot 2: PAM50
# -----------------------------
p_pam50 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = PAM50)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_color_manual(values = COLORS$PAM50, na.value = "grey80") +
  labs(
    title = "PCA of TCGA-BRCA transcriptomes",
    subtitle = "Top 5,000 variable genes • Colored by PAM50 subtype",
    x = pc1_lab,
    y = pc2_lab,
    color = "PAM50"
  ) +
  theme_project(BASE_SIZE)

out_base_pam50 <- file.path(out_dir, "PCA_expr1218_PAM50")
save_figure(p_pam50, out_base_pam50)

# -----------------------------
# Write QC report
# -----------------------------
qc_file <- file.path(out_dir, "PCA_QC_summary_expr1218.txt")
writeLines(qc_report, qc_file)

# -----------------------------
# Final messages
# -----------------------------
cat("\n🏁 STEP 4 COMPLETED SUCCESSFULLY\n")
cat("Figures saved in:\n", out_dir, "\n")
cat("Variance table:", var_file, "\n")
cat("QC report:", qc_file, "\n")
cat("Outliers:", outlier_file, "\n")
