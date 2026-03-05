# ================================
# STEP 11 — FINAL TARGET INTERSECTION
# ================================

cat("🧬 Starting FINAL TARGET INTERSECTION...\n")

base_dir <- "D:/TNBC_Project_CLEAN/results"

# ------------------------
# 1. Load ML ranking
# ------------------------
ml_file <- file.path(base_dir, "ML_Target_Ranking", "Final_target_ranking.tsv")
stopifnot(file.exists(ml_file))

ml <- read.delim(ml_file, stringsAsFactors = FALSE)
ml_genes <- unique(ml$Gene)

cat("✅ ML genes:", length(ml_genes), "\n")

# ------------------------
# 2. Load STRING hub genes
# ------------------------
string_file <- file.path(base_dir, "STRING_hub_genes.txt")
stopifnot(file.exists(string_file))

string_genes <- readLines(string_file)
string_genes <- unique(trimws(string_genes))

cat("✅ STRING hub genes:", length(string_genes), "\n")

# ------------------------
# 3. Load DE genes
# ------------------------

de_file <- "D:/TNBC_Project_CLEAN/results/DE_TNBC_vs_nonTNBC/DE_full_TNBC_vs_nonTNBC.tsv"

stopifnot(file.exists(de_file))

de <- read.delim(de_file, stringsAsFactors = FALSE)

# use adj.P.Val and logFC
de_sig <- subset(de, adj.P.Val < 0.05 & abs(logFC) > 1)
de_genes <- unique(de_sig$Gene)

cat("✅ DE genes:", length(de_genes), "\n")

# ------------------------
# 4. Load Mutation genes
# ------------------------
mut_file <- "D:/TNBC_Project_CLEAN/results/MUT_TNBC_vs_nonTNBC/Mutation_comparison_significant.tsv"

stopifnot(file.exists(mut_file))

mut <- read.delim(mut_file, stringsAsFactors = FALSE)
mut_genes <- unique(mut$Gene)

cat("✅ Mutated genes:", length(mut_genes), "\n")

# ------------------------
# 5. Union of all candidates
# ------------------------
all_genes <- unique(c(ml_genes, de_genes, mut_genes, string_genes))

# ------------------------
# 6. Build evidence table
# ------------------------
evidence <- data.frame(
  Gene = all_genes,
  In_ML = all_genes %in% ml_genes,
  In_DE = all_genes %in% de_genes,
  In_MUT = all_genes %in% mut_genes,
  In_STRING = all_genes %in% string_genes
)

evidence$Support_Count <- rowSums(evidence[,2:5])

# ------------------------
# 7. Select high-confidence targets
# ------------------------
final_targets <- subset(evidence, Support_Count >= 3)

# ------------------------
# 8. Save results
# ------------------------
outdir <- file.path(base_dir, "Final_Targets")
dir.create(outdir, showWarnings = FALSE)

write.table(evidence,
            file.path(outdir, "All_candidates_with_evidence.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(final_targets,
            file.path(outdir, "Final_high_confidence_targets.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------
# 9. Print summary
# ------------------------
cat("\n==============================\n")
cat("🎯 FINAL TARGET SELECTION SUMMARY\n")
cat("==============================\n")
cat("Total candidate genes:", nrow(evidence), "\n")
cat("High-confidence targets (>=3 evidences):", nrow(final_targets), "\n\n")

print(final_targets[order(-final_targets$Support_Count), ])

cat("\n✅ Final targets saved to:", outdir, "\n")
cat("🏁 STEP 11 COMPLETE!\n")
