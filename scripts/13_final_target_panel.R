############################################
# STEP 10 — Final Target Panel & Drug Mapping
############################################

cat("🧬 Starting STEP 10: Final therapeutic target panel...\n")

library(data.table)

base_dir <- "D:/TNBC_Project_CLEAN"

ml_file   <- file.path(base_dir, "results/ML_Target_Ranking/Final_target_ranking.tsv")
hub_file  <- file.path(base_dir, "results/STRING_Network/Top20_STRING_hubs.tsv")
mut_file  <- file.path(base_dir, "results/MUT_TNBC_vs_nonTNBC/Mutation_comparison_significant.tsv")
de_file   <- file.path(base_dir, "results/DE_TNBC_vs_nonTNBC/DE_full_TNBC_vs_nonTNBC.tsv")

out_dir <- file.path(base_dir, "results/Final_Targets")
dir.create(out_dir, showWarnings = FALSE)

# ----------------------------
# Load data
# ----------------------------
ml  <- fread(ml_file)
hub <- fread(hub_file)
mut <- fread(mut_file)
de  <- fread(de_file)

# ----------------------------
# Select strong DE genes
# ----------------------------
de_sig <- de[adj.P.Val < 0.01 & abs(logFC) > 1.5]

# ----------------------------
# Build evidence table
# ----------------------------
genes <- unique(c(
  head(ml$Gene, 50),
  hub$Gene,
  mut$Gene,
  head(de_sig$Gene, 200)
))

final <- data.table(Gene = genes)

final[, ML := Gene %in% ml$Gene[1:50]]
final[, Hub := Gene %in% hub$Gene]
final[, Mutation := Gene %in% mut$Gene]
final[, DE := Gene %in% de_sig$Gene]

final[, EvidenceScore := ML + Hub + Mutation + DE]

# ----------------------------
# Keep strong candidates
# ----------------------------
final <- final[EvidenceScore >= 2]
final <- final[order(-EvidenceScore)]

# ----------------------------
# Add known drugs manually (curated)
# ----------------------------
drug_map <- data.table(
  Gene = c("AURKA","PLK1","CDK1","TOP2A","EGFR","KIF11","TTK","ERBB2"),
  Drug = c("Alisertib","Volasertib","RO-3306","Doxorubicin","Erlotinib","Ispinesib","BAY-1217389","Lapatinib")
)

final <- merge(final, drug_map, by = "Gene", all.x = TRUE)

# ----------------------------
# Save
# ----------------------------
fwrite(final, file.path(out_dir, "Final_Therapeutic_Target_Panel.tsv"), sep="\t")

cat("✅ Final targets:", nrow(final), "\n")
cat("📁 Saved to:", out_dir, "\n")
