############################################################
# QC SAMPLE FLOW SUMMARY SCRIPT (CLEAN & FINAL)
# Purpose:
#   Show how 1231 samples became 1218 samples after QC
############################################################

cat("\n========== QC SAMPLE FLOW REPORT ==========\n\n")

library(SummarizedExperiment)

# ---- Files ----
raw_file   <- "D:/TNBC_Project_CLEAN/data_raw/TCGA_BRCA_STAR_counts_FINAL.rds"
clean_file <- "D:/TNBC_Project_CLEAN/data_processed/expr_clean.rds"

stopifnot(file.exists(raw_file))
stopifnot(file.exists(clean_file))

# ---- Load datasets ----
cat("Loading RAW dataset...\n")
se_raw <- readRDS(raw_file)

cat("Loading CLEAN dataset...\n")
expr_clean <- readRDS(clean_file)

# ---- Get dimensions ----
raw_genes   <- nrow(assay(se_raw))
raw_samples <- ncol(assay(se_raw))

clean_genes   <- nrow(expr_clean)
clean_samples <- ncol(expr_clean)

removed_samples <- raw_samples - clean_samples

# ---- Print report ----
cat("RAW DATASET:\n")
cat("  Genes   :", raw_genes, "\n")
cat("  Samples :", raw_samples, "\n\n")

cat("CLEAN DATASET:\n")
cat("  Genes   :", clean_genes, "\n")
cat("  Samples :", clean_samples, "\n\n")

cat("SAMPLE FLOW SUMMARY:\n")
cat("  Downloaded samples :", raw_samples, "\n")
cat("  Final samples      :", clean_samples, "\n")
cat("  Removed samples    :", removed_samples, "\n\n")

# ---- Save summary table ----
qc_summary <- data.frame(
  Metric = c("Downloaded samples", "Final samples", "Removed samples"),
  Count  = c(raw_samples, clean_samples, removed_samples)
)

out_file <- "D:/TNBC_Project_CLEAN/results/QC_sample_flow_summary.tsv"
write.table(qc_summary, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("QC summary saved to:\n", out_file, "\n")
cat("\n========== QC REPORT COMPLETED ==========\n")

h
