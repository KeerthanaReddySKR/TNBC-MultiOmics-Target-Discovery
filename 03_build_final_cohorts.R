# =========================================
# 03_build_final_cohorts.R
# =========================================

rm(list = ls())
cat("🧬 Building final frozen cohorts...\n")

# -------- Paths --------
expr_file   <- "D:/TNBC_Project_CLEAN/data_processed/expr_clean.rds"
labels_file <- "D:/TNBC_Project_CLEAN/data_processed/BRCA_TNBC_labels.rds"
clin_file   <- "D:/TNBC_Project_CLEAN/data_raw/TCGA_BRCA_clinical_Xena.tsv"
maf_file    <- "D:/TNBC_Project_CLEAN/data_raw/TCGA_BRCA_mutation_maf_object.rds"

out_1218 <- "D:/TNBC_Project_CLEAN/data_processed/expr_1218.rds"
out_1024 <- "D:/TNBC_Project_CLEAN/data_processed/expr_1024.rds"

# -------- Load libraries --------
suppressPackageStartupMessages({
  library(data.table)
})

# -------- Load data --------
cat("📥 Loading expression...\n")
expr <- readRDS(expr_file)

cat("📥 Loading labels...\n")
labels <- readRDS(labels_file)

cat("📥 Loading clinical...\n")
clin <- fread(clin_file)

cat("📥 Loading mutation...\n")
maf <- readRDS(maf_file)

# -------- Align expression + labels + clinical --------
common1 <- Reduce(intersect, list(
  colnames(expr),
  labels$sampleID,
  clin$sampleID
))

cat("Samples common to expr + labels + clinical:", length(common1), "\n")

# Subset
expr1 <- expr[, common1]
labels1 <- labels[match(common1, sampleID)]
clin1 <- clin[match(common1, sampleID)]

# -------- Remove samples with NA TNBC label --------
keep <- !is.na(labels1$is_TNBC)

cat("Removing samples with NA TNBC label:", sum(!keep), "\n")

expr1 <- expr1[, keep]
labels1 <- labels1[keep]
clin1 <- clin1[keep]

cat("Final samples for expr_1218:", ncol(expr1), "\n")
cat("Genes:", nrow(expr1), "\n")

# -------- Save expr_1218 --------
expr_1218 <- list(
  expr = expr1,
  labels = labels1,
  clinical = clin1
)

saveRDS(expr_1218, out_1218)
cat("💾 Saved expr_1218 to:", out_1218, "\n")

# -------- Build expr_1024 (with mutation) --------

# Get mutation sample IDs (first 15 chars)
maf_samples <- substr(maf@clinical.data$Tumor_Sample_Barcode, 1, 15)
maf_samples <- unique(maf_samples)

common2 <- intersect(colnames(expr1), maf_samples)

cat("Samples common to expr_1218 + mutation:", length(common2), "\n")

# Subset again
expr2 <- expr1[, common2]
labels2 <- labels1[match(common2, labels1$sampleID)]
clin2 <- clin1[match(common2, clin1$sampleID)]

cat("Final samples for expr_1024:", ncol(expr2), "\n")
cat("Genes:", nrow(expr2), "\n")

# -------- Save expr_1024 --------
expr_1024 <- list(
  expr = expr2,
  labels = labels2,
  clinical = clin2,
  maf = maf
)

saveRDS(expr_1024, out_1024)
cat("💾 Saved expr_1024 to:", out_1024, "\n")

# -------- Final report --------
cat("\n🏁 FINAL FROZEN COHORTS SUMMARY\n")
cat("expr_1218: genes =", nrow(expr1), " samples =", ncol(expr1), "\n")
cat("expr_1024: genes =", nrow(expr2), " samples =", ncol(expr2), "\n")

cat("\n🎉 Cohort construction completed successfully.\n")
