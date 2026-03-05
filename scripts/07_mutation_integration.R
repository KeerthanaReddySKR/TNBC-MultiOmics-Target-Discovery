# ============================================================
# STEP 7: Somatic mutation integration (TNBC vs non-TNBC)
# FINAL BULLETPROOF VERSION (NO mafCompare)
# ============================================================

rm(list = ls())
cat("🧬 Starting STEP 7: Mutation integration...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
in_file <- "D:/TNBC_Project_CLEAN/data_processed/expr_1024.rds"
out_dir <- "D:/TNBC_Project_CLEAN/results/MUT_TNBC_vs_nonTNBC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load integrated object
# -----------------------------
cat("📥 Loading expr_1024 object...\n")
obj <- readRDS(in_file)

expr     <- obj$expr
labels   <- obj$labels
clinical <- obj$clinical
maf      <- obj$maf

cat("Expression samples:", ncol(expr), "\n")
cat("MAF samples (raw full barcodes):", length(unique(maf@data$Tumor_Sample_Barcode)), "\n")

stopifnot("is_TNBC" %in% colnames(labels))

# -----------------------------
# Build 12-char TCGA IDs from expression
# -----------------------------
labels$SampleShort <- substr(labels$sampleID, 1, 12)

# -----------------------------
# Harmonize MAF barcodes to 12-char patient IDs
# -----------------------------
cat("🔧 Harmonizing MAF barcodes to 12-char patient IDs...\n")

maf@data$Tumor_Sample_Barcode_full <- maf@data$Tumor_Sample_Barcode
maf@data$Tumor_Sample_Barcode <- substr(maf@data$Tumor_Sample_Barcode, 1, 12)

maf_ids_short <- unique(maf@data$Tumor_Sample_Barcode)
cat("Unique MAF patient IDs after harmonization:", length(maf_ids_short), "\n")

# -----------------------------
# Build TNBC / non-TNBC lists
# -----------------------------
tnbc_samples    <- labels$SampleShort[labels$is_TNBC == TRUE]
nontnbc_samples <- labels$SampleShort[labels$is_TNBC == FALSE]

cat("TNBC samples (expression):", length(tnbc_samples), "\n")
cat("Non-TNBC samples (expression):", length(nontnbc_samples), "\n")

# -----------------------------
# Overlap check
# -----------------------------
overlap_tnbc    <- length(intersect(tnbc_samples, maf_ids_short))
overlap_nontnbc <- length(intersect(nontnbc_samples, maf_ids_short))

cat("Overlap TNBC vs MAF:", overlap_tnbc, "\n")
cat("Overlap non-TNBC vs MAF:", overlap_nontnbc, "\n")

if (overlap_tnbc == 0 || overlap_nontnbc == 0) {
  stop("❌ FATAL: No overlap between expression cohort and MAF after harmonization.")
}

# -----------------------------
# Subset MAF
# -----------------------------
cat("🔪 Subsetting MAF into TNBC and non-TNBC...\n")

maf_tnbc <- subsetMaf(maf = maf, tsb = tnbc_samples, mafObj = TRUE)
maf_nontnbc <- subsetMaf(maf = maf, tsb = nontnbc_samples, mafObj = TRUE)

n_TNBC <- length(unique(maf_tnbc@data$Tumor_Sample_Barcode))
n_nonTNBC <- length(unique(maf_nontnbc@data$Tumor_Sample_Barcode))

cat("MAF TNBC samples with mutations:", n_TNBC, "\n")
cat("MAF non-TNBC samples with mutations:", n_nonTNBC, "\n")

# -----------------------------
# Oncoplots (Top 20)
# -----------------------------
cat("🎨 Generating oncoplots...\n")

pdf(file.path(out_dir, "Oncoplot_TNBC_Top20.pdf"), width = 10, height = 8)
oncoplot(maf = maf_tnbc, top = 20, removeNonMutated = TRUE, titleText = "TNBC: Top 20 Mutated Genes")
dev.off()

png(file.path(out_dir, "Oncoplot_TNBC_Top20.png"), width = 3000, height = 2400, res = 300)
oncoplot(maf = maf_tnbc, top = 20, removeNonMutated = TRUE, titleText = "TNBC: Top 20 Mutated Genes")
dev.off()

pdf(file.path(out_dir, "Oncoplot_nonTNBC_Top20.pdf"), width = 10, height = 8)
oncoplot(maf = maf_nontnbc, top = 20, removeNonMutated = TRUE, titleText = "Non-TNBC: Top 20 Mutated Genes")
dev.off()

png(file.path(out_dir, "Oncoplot_nonTNBC_Top20.png"), width = 3000, height = 2400, res = 300)
oncoplot(maf = maf_nontnbc, top = 20, removeNonMutated = TRUE, titleText = "Non-TNBC: Top 20 Mutated Genes")
dev.off()

# -----------------------------
# MANUAL, BULLETPROOF DIFFERENTIAL MUTATION TEST
# -----------------------------
cat("🔬 Building mutation count table manually...\n")

tnbc_mut <- table(maf_tnbc@data$Hugo_Symbol)
nontnbc_mut <- table(maf_nontnbc@data$Hugo_Symbol)

all_genes <- sort(unique(c(names(tnbc_mut), names(nontnbc_mut))))

mut_table <- data.frame(
  Gene = all_genes,
  TNBC_mut = as.integer(tnbc_mut[all_genes]),
  NonTNBC_mut = as.integer(nontnbc_mut[all_genes])
)

mut_table$TNBC_mut[is.na(mut_table$TNBC_mut)] <- 0
mut_table$NonTNBC_mut[is.na(mut_table$NonTNBC_mut)] <- 0

# Add wildtype counts
mut_table$TNBC_wt <- n_TNBC - mut_table$TNBC_mut
mut_table$NonTNBC_wt <- n_nonTNBC - mut_table$NonTNBC_mut

# Filter untestable genes
mut_table <- mut_table[(mut_table$TNBC_mut + mut_table$NonTNBC_mut) >= 5, ]

cat("Genes tested after filtering:", nrow(mut_table), "\n")

# Run Fisher safely
cat("📊 Running Fisher tests safely...\n")

fisher_results <- apply(mut_table, 1, function(x) {
  mat <- matrix(
    c(
      as.numeric(x["TNBC_mut"]),
      as.numeric(x["TNBC_wt"]),
      as.numeric(x["NonTNBC_mut"]),
      as.numeric(x["NonTNBC_wt"])
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  if (any(is.na(mat)) || any(mat < 0)) {
    return(c(p = NA, OR = NA))
  }
  
  ft <- fisher.test(mat)
  c(p = ft$p.value, OR = as.numeric(ft$estimate))
})

fisher_results <- t(fisher_results)

mut_table$p_value <- fisher_results[, "p"]
mut_table$OR <- fisher_results[, "OR"]
mut_table$FDR <- p.adjust(mut_table$p_value, method = "fdr")
mut_table$log2OR <- log2(mut_table$OR)

sig_table <- mut_table[!is.na(mut_table$FDR) & mut_table$FDR < 0.05, ]

cat("Significant differential mutated genes (FDR < 0.05):", nrow(sig_table), "\n")

# Save tables
fwrite(mut_table, file.path(out_dir, "Mutation_comparison_full.tsv"), sep = "\t")
fwrite(sig_table, file.path(out_dir, "Mutation_comparison_significant.tsv"), sep = "\t")

# -----------------------------
# Barplot of top hits
# -----------------------------
cat("🎨 Plotting differential mutation barplot...\n")

if (nrow(sig_table) > 0) {
  
  top_plot <- sig_table[order(sig_table$FDR), ]
  top_plot <- top_plot[1:min(10, nrow(top_plot)), ]
  
  p <- ggplot(top_plot, aes(x = reorder(Gene, log2OR), y = log2OR, fill = log2OR > 0)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = "none") +
    labs(
      title = "Differentially mutated genes (TNBC vs non-TNBC)",
      subtitle = "log2 Odds Ratio (Fisher test, FDR < 0.05)",
      x = "",
      y = "log2 Odds Ratio (TNBC / non-TNBC)"
    ) +
    theme_classic(base_size = 15)
  
  ggsave(file.path(out_dir, "Differential_Mutation_Barplot.pdf"), p, width = 7, height = 6)
  ggsave(file.path(out_dir, "Differential_Mutation_Barplot.png"), p, width = 7, height = 6, dpi = 300)
}

# -----------------------------
# QC Summary
# -----------------------------
qc_file <- file.path(out_dir, "MUT_QC_summary.txt")

qc_lines <- c(
  "STEP 7: Mutation Integration QC Summary",
  "======================================",
  paste("Date:", Sys.time()),
  "",
  paste("Total samples in expr_1024:", ncol(expr)),
  paste("TNBC samples (expression):", length(tnbc_samples)),
  paste("Non-TNBC samples (expression):", length(nontnbc_samples)),
  "",
  paste("MAF TNBC samples with mutations:", n_TNBC),
  paste("MAF non-TNBC samples with mutations:", n_nonTNBC),
  "",
  paste("Genes tested (after >=5 total mutations filter):", nrow(mut_table)),
  paste("Significant differential mutated genes (FDR < 0.05):", nrow(sig_table)),
  "",
  "Top significant genes:",
  if (nrow(sig_table) > 0) paste0("  - ", sig_table$Gene[1:min(10, nrow(sig_table))]) else "  (none)"
)

writeLines(qc_lines, qc_file)

# -----------------------------
# Finish
# -----------------------------
cat("✅ STEP 7 Mutation integration completed successfully.\n")
cat("Results written to:", out_dir, "\n")
