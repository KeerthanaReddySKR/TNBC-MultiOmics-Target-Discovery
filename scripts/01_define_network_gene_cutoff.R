# ================================
# STEP 1: Define objective cutoff for STRING network
# ================================

library(ggplot2)

in_file <- "D:/TNBC_Project_CLEAN/results/ML_Target_Ranking/Final_target_ranking.tsv"
out_dir <- "D:/TNBC_Project_CLEAN/results/ML_Target_Ranking/"

stopifnot(file.exists(in_file))

df <- read.delim(in_file, stringsAsFactors = FALSE)

# Sort by score
df <- df[order(-df$Score), ]

# ---- NORMALIZE SCORE (CRITICAL FIX) ----
score_min <- min(df$Score)
score_max <- max(df$Score)

df$Score_norm <- (df$Score - score_min) / (score_max - score_min)

# Replace zeros if any
df$Score_norm[df$Score_norm < 1e-6] <- 1e-6

# Cumulative sum
df$cum_norm <- cumsum(df$Score_norm) / sum(df$Score_norm)

# Find cutoffs
cut85 <- which(df$cum_norm >= 0.85)[1]
cut90 <- which(df$cum_norm >= 0.90)[1]

cat("Genes needed for 85% signal:", cut85, "\n")
cat("Genes needed for 90% signal:", cut90, "\n")

# Save gene lists
write.table(
  df[1:cut85, c("Gene", "Score", "Score_norm", "cum_norm")],
  file = paste0(out_dir, "Genes_for_STRING_85pct.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  df[1:cut90, c("Gene", "Score", "Score_norm", "cum_norm")],
  file = paste0(out_dir, "Genes_for_STRING_90pct.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Plot curve
p <- ggplot(df, aes(x = seq_along(cum_norm), y = cum_norm)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = cut85, linetype = "dotted", color = "red") +
  geom_vline(xintercept = cut90, linetype = "dotted", color = "darkgreen") +
  labs(
    title = "Cumulative contribution of ranked genes (normalized score)",
    x = "Top N genes",
    y = "Cumulative normalized score"
  ) +
  theme_minimal()

ggsave(
  paste0(out_dir, "Cumulative_score_curve.png"),
  p,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved plots and gene lists to:", out_dir, "\n")

