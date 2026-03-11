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

tnbc_file  <- "D:/TNBC_Project_CLEAN/data_processed/TNBC_123_patients.tsv"

# -------- Libraries --------
suppressPackageStartupMessages({
  library(data.table)
})

# -------- Load expression --------
cat("📥 Loading expression...\n")
expr <- readRDS(expr_file)

# -------- Load clinical --------
cat("📥 Loading clinical...\n")
clin <- fread(clin_file)

# -------- Match samples --------
common <- intersect(colnames(expr), clin$sampleID)

cat("Common samples (expr ∩ clinical):", length(common), "\n")

if (length(common) < 900) {
  stop("❌ ERROR: Too few overlapping samples. Something is wrong.")
}

clin2 <- clin[match(common, sampleID)]

# -------- Extract receptors --------
ER   <- clin2$ER_Status_nature2012
PR   <- clin2$PR_Status_nature2012
HER2 <- clin2$HER2_Final_Status_nature2012

# -------- Normalize receptor text --------
ER   <- toupper(trimws(ER))
PR   <- toupper(trimws(PR))
HER2 <- toupper(trimws(HER2))

# -------- Define negative receptors --------
is_ER_neg   <- ER   %in% c("NEGATIVE", "NEG")
is_PR_neg   <- PR   %in% c("NEGATIVE", "NEG")
is_HER2_neg <- HER2 %in% c("NEGATIVE", "NEG")

# -------- Define TNBC --------
is_TNBC <- is_ER_neg & is_PR_neg & is_HER2_neg

# -------- Handle missing values --------
valid <- !is.na(is_ER_neg) & !is.na(is_PR_neg) & !is.na(is_HER2_neg)

cat("Samples with complete receptor info:", sum(valid), "\n")

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

# Add readable class column
labels[, class := ifelse(is_TNBC == TRUE, "TNBC", "NON_TNBC")]

# -------- Sanity check --------
cat("\n🔎 Verifying TNBC receptor logic...\n")

tnbc_check <- labels[is_TNBC == TRUE, .(ER, PR, HER2)]

print(unique(tnbc_check))

if (!all(tnbc_check$ER == "NEGATIVE" &
         tnbc_check$PR == "NEGATIVE" &
         tnbc_check$HER2 == "NEGATIVE")) {
  
  stop("❌ ERROR: Some TNBC samples are not triple negative!")
  
} else {
  
  cat("✅ All TNBC samples confirmed ER-, PR-, HER2-\n")
}

# -------- Show example patients --------
cat("\n🧪 Example TNBC patients:\n")
print(labels[class == "TNBC"][1:10])

cat("\n🧪 Example non-TNBC patients:\n")
print(labels[class == "NON_TNBC"][1:10])

# -------- Extract 123 TNBC patients --------
tnbc_patients <- labels[class == "TNBC", .(sampleID, ER, PR, HER2)]

cat("\nTotal TNBC patients:", nrow(tnbc_patients), "\n")

# -------- Save TNBC patient list --------
fwrite(tnbc_patients, tnbc_file, sep="\t")

cat("💾 Saved TNBC patient list to:\n", tnbc_file, "\n")

# -------- Save full label table --------
saveRDS(labels, out_labels)
fwrite(labels, out_table, sep = "\t")

cat("\n💾 Saved full labels to:\n", out_labels, "\n", out_table, "\n")

# -------- TNBC proportion --------
tnbc_rate <- sum(labels$is_TNBC, na.rm=TRUE) / sum(valid)

cat("\n📊 TNBC proportion in dataset:",
    round(tnbc_rate*100,2), "%\n")

cat("\n🎉 TNBC label construction completed successfully.\n")