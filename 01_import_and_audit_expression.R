# ================================
# 01_import_and_audit_expression.R
# ================================

rm(list = ls())
cat("🧬 Starting expression import & audit...\n")

# -------- Paths --------
expr_file <- "D:/TNBC_Project_CLEAN/data_raw/TCGA.BRCA.sampleMap_HiSeqV2.gz"
clin_file <- "D:/TNBC_Project_CLEAN/data_raw/TCGA_BRCA_clinical_Xena.tsv"

out_file  <- "D:/TNBC_Project_CLEAN/data_processed/expr_clean.rds"

# -------- Load libraries --------
suppressPackageStartupMessages({
  library(data.table)
})

# -------- Load expression --------
cat("📥 Loading expression matrix...\n")
expr_dt <- fread(expr_file)

cat("Dimensions of raw table:", dim(expr_dt), "\n")

# First column is gene ID
gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids

# -------- Basic checks --------
cat("Genes:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n")

cat("\n🧪 First 5 sample IDs:\n")
print(colnames(expr_mat)[1:5])

# Check TCGA format
cat("\n🧪 Do sample IDs look like TCGA barcodes?\n")
print(substr(colnames(expr_mat)[1:5], 1, 12))

# -------- Check duplicates --------
cat("\nDuplicate genes:", sum(duplicated(rownames(expr_mat))), "\n")
cat("Duplicate samples:", sum(duplicated(colnames(expr_mat))), "\n")

# -------- Load clinical --------
cat("\n📥 Loading clinical data...\n")
clin <- fread(clin_file)

cat("Clinical dimensions:", dim(clin), "\n")
cat("First 5 clinical sampleIDs:\n")
print(head(clin$sampleID, 5))

# -------- Intersection check --------
common <- intersect(colnames(expr_mat), clin$sampleID)

cat("\n🧪 Number of common samples between expression and clinical:\n")
cat(length(common), "\n")

# -------- Verdict --------
if (length(common) < 900) {
  stop("❌ ERROR: Too few overlapping samples. Something is wrong.")
} else {
  cat("✅ PASS: Expression and clinical are properly linkable.\n")
}

# -------- Subset to common samples --------
expr_mat2 <- expr_mat[, common]

cat("\nAfter subsetting:\n")
cat("Genes:", nrow(expr_mat2), "\n")
cat("Samples:", ncol(expr_mat2), "\n")

# -------- Save clean object --------
saveRDS(expr_mat2, out_file)
cat("\n💾 Saved clean expression to:\n", out_file, "\n")

cat("\n🎉 Expression import & audit completed successfully.\n")
