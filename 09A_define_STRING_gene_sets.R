# ============================================================
# STEP 9A — Define STRING input gene sets from ML ranking
# ============================================================

cat("🧬 Defining STRING input gene sets from ML ranking...\n")

base_dir <- "D:/TNBC_Project_CLEAN"
ml_file <- file.path(base_dir, "results/ML_Target_Ranking/Final_target_ranking.tsv")

stopifnot(file.exists(ml_file))

out_dir <- file.path(base_dir, "results/STRING_Input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load ML ranking
ml <- read.table(ml_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Sanity check
stopifnot("Gene" %in% colnames(ml))

# Sort by final rank (best first)
if ("Final_rank" %in% colnames(ml)) {
  ml <- ml[order(ml$Final_rank), ]
}

# Define cutoffs
cuts <- c(100, 200, 300)

for (k in cuts) {
  genes <- unique(ml$Gene[1:k])
  
  out_file <- file.path(out_dir, paste0("STRING_Top", k, ".txt"))
  writeLines(genes, out_file)
  
  cat("✅ Saved:", out_file, " (n =", length(genes), ")\n")
}

cat("🎯 DONE. Now upload these to STRING one by one and check network quality.\n")
