# =========================================
# 02_define_TNBC_and_build_labels.R
# =========================================

rm(list = ls())
cat("🧬 Defining TNBC vs non-TNBC labels...\n")

# -------- Paths --------
expr_file <- "D:/TNBC_Project_CLEAN/data_processed/expr_clean.rds"
clin_file <- "D:/TNBC_Project_CLEAN/data_raw/TCGA_BRCA_clinical_Xena.tsv"

out_labels <- "D:/TNBC_Project_CLEAN/data_processed/BRCA_TNBC_labels.rds"
out_table  <- "D:/TNBC_Project_CLEAN/data_processed/BRCA_TNBC_labels.tsv"

# -------- Load libraries --------
suppressPackageStartupMessages({
  library(data.table)
})

# -------- Load data --------
cat("📥 Loading expression...\n")
expr <- readRDS(expr_file)

cat("📥 Loading clinical...\n")
clin <- fread(clin_file)

# -------- Keep only samples present in expression --------
common <- intersect(colnames(expr), clin$sampleID)

cat("Common samples (expr ∩ clinical):", length(common), "\n")

if (length(common) < 900) {
  stop("❌ Too few overlapping samples. Something is wrong.")
}

clin2 <- clin[match(common, sampleID)]

# -------- Extract receptor status --------
ER   <- clin2$ER_Status_nature2012
PR   <- clin2$PR_Status_nature2012
HER2 <- clin2$HER2_Final_Status_nature2012

# -------- Normalize text (safety) --------
ER   <- toupper(trimws(ER))
PR   <- toupper(trimws(PR))
HER2 <- toupper(trimws(HER2))

# -------- Define negativity --------
is_ER_neg   <- ER   %in% c("NEGATIVE", "NEG")
is_PR_neg   <- PR   %in% c("NEGATIVE", "NEG")
is_HER2_neg <- HER2 %in% c("NEGATIVE", "NEG")

# -------- Define TNBC --------
is_TNBC <- is_ER_neg & is_PR_neg & is_HER2_neg

# -------- Handle missing values --------
valid <- !is.na(is_ER_neg) & !is.na(is_PR_neg) & !is.na(is_HER2_neg)

cat("Samples with complete receptor info:", sum(valid), "\n")

# Set TNBC to NA if any receptor is missing
is_TNBC[!valid] <- NA

# -------- Summary --------
cat("\n🧪 TNBC label summary (including NAs):\n")
print(table(is_TNBC, useNA = "ifany"))

cat("\n🧪 TNBC label summary (complete cases only):\n")
print(table(is_TNBC[valid]))

# -------- Build label table --------
labels <- data.table(
  sampleID = clin2$sampleID,
  ER = ER,
  PR = PR,
  HER2 = HER2,
  is_TNBC = is_TNBC
)

# -------- Save --------
saveRDS(labels, out_labels)
fwrite(labels, out_table, sep = "\t")

cat("\n💾 Saved labels to:\n", out_labels, "\n", out_table, "\n")

cat("\n🎉 TNBC label construction completed successfully.\n")
