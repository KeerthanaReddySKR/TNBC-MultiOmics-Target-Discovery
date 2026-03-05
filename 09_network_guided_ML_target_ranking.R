# ============================================================
# STEP 9: Network-guided Machine Learning & Target Ranking
# ============================================================

rm(list = ls())
cat("🧬 Starting STEP 9: Network-guided ML & target ranking...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(glmnet)
  library(pROC)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
base_dir <- "D:/TNBC_Project_CLEAN"

expr_file <- file.path(base_dir, "data_processed/expr_1218.rds")

de_dir   <- file.path(base_dir, "results/DE_TNBC_vs_nonTNBC")
net_dir  <- file.path(base_dir, "results/STRING_Network")
mut_dir  <- file.path(base_dir, "results/MUT_TNBC_vs_nonTNBC")

out_dir <- file.path(base_dir, "results/ML_Target_Ranking")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(expr_file), dir.exists(de_dir), dir.exists(net_dir), dir.exists(mut_dir))

# -----------------------------
# Auto-detect files
# -----------------------------
cat("🔎 Auto-detecting input files...\n")

de_file  <- list.files(de_dir, pattern = "DE_.*full.*\\.tsv$", full.names = TRUE)[1]
hub_file <- list.files(net_dir, pattern = "Top.*hub.*\\.tsv$", full.names = TRUE)[1]
mut_file <- list.files(mut_dir, pattern = "Mutation_.*significant.*\\.tsv$", full.names = TRUE)[1]

if (is.na(de_file) || is.na(hub_file) || is.na(mut_file)) {
  stop("❌ Could not auto-detect one or more input files. Please check result folders.")
}

cat("✅ Using files:\n")
cat("DE:   ", de_file, "\n")
cat("Hubs: ", hub_file, "\n")
cat("Mut:  ", mut_file, "\n")

# -----------------------------
# Load data
# -----------------------------
cat("📥 Loading expr_1218...\n")
obj <- readRDS(expr_file)
expr <- obj$expr
labels <- obj$labels

y <- as.factor(ifelse(labels$is_TNBC, "TNBC", "nonTNBC"))
names(y) <- labels$sampleID

# -----------------------------
# Load feature sources
# -----------------------------
cat("📥 Loading DE, network hubs, mutation genes...\n")

de  <- fread(de_file)
hub <- fread(hub_file)
mut <- fread(mut_file)

# -----------------------------
# Standardize column names
# -----------------------------
if (!"Gene" %in% colnames(de))  setnames(de, 1, "Gene")
if (!"Gene" %in% colnames(hub)) setnames(hub, 1, "Gene")
if (!"Gene" %in% colnames(mut)) setnames(mut, 1, "Gene")

# -----------------------------
# Select feature genes
# -----------------------------

# Top DE genes
de_sig <- de[adj.P.Val < 0.01 & abs(logFC) > 1.5]
de_genes <- de_sig$Gene

# Network hubs
hub_genes <- hub$Gene

# Mutated genes
mut_genes <- mut$Gene

# Union
feature_genes <- unique(c(de_genes, hub_genes, mut_genes))

cat("DE genes used:", length(de_genes), "\n")
cat("Hub genes used:", length(hub_genes), "\n")
cat("Mut genes used:", length(mut_genes), "\n")
cat("Total unique features:", length(feature_genes), "\n")

# Keep only genes present in expression matrix
feature_genes <- intersect(feature_genes, rownames(expr))
cat("Features present in expression matrix:", length(feature_genes), "\n")

# -----------------------------
# Build ML matrix
# -----------------------------
cat("🧱 Building ML feature matrix...\n")

X <- t(expr[feature_genes, ])
X <- scale(X)

# Ensure sample order matches labels
X <- X[names(y), ]
stopifnot(all(rownames(X) == names(y)))

# -----------------------------
# Random Forest
# -----------------------------
cat("🌲 Training Random Forest...\n")

set.seed(123)
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)

rf_prob <- predict(rf, X, type = "prob")[, "TNBC"]

roc_rf <- roc(response = y, predictor = rf_prob, levels = c("nonTNBC", "TNBC"))
auc_rf <- auc(roc_rf)

cat("RF AUC:", auc_rf, "\n")

# RF importance
rf_imp <- importance(rf, type = 2)
rf_imp_df <- data.frame(
  Gene = rownames(rf_imp),
  RF_Importance = rf_imp[, 1]
)
rf_imp_df <- rf_imp_df[order(-rf_imp_df$RF_Importance), ]

fwrite(rf_imp_df, file.path(out_dir, "Feature_importance_RF.tsv"), sep = "\t")

# -----------------------------
# LASSO (glmnet)
# -----------------------------
cat("📐 Training LASSO (glmnet)...\n")

y_bin <- ifelse(y == "TNBC", 1, 0)

cvfit <- cv.glmnet(as.matrix(X), y_bin, family = "binomial", alpha = 1)
best_lambda <- cvfit$lambda.min

cat("Best lambda:", best_lambda, "\n")

lasso <- glmnet(as.matrix(X), y_bin, family = "binomial", alpha = 1, lambda = best_lambda)

coef_mat <- coef(lasso)
coef_df <- data.frame(
  Gene = rownames(coef_mat),
  LASSO_Coefficient = as.numeric(coef_mat)
)

coef_df <- coef_df[coef_df$Gene != "(Intercept)", ]
coef_df <- coef_df[coef_df$LASSO_Coefficient != 0, ]
coef_df <- coef_df[order(-abs(coef_df$LASSO_Coefficient)), ]

fwrite(coef_df, file.path(out_dir, "Feature_importance_LASSO.tsv"), sep = "\t")

# -----------------------------
# ROC plot
# -----------------------------
cat("📈 Saving ROC plot...\n")

roc_df <- data.frame(
  FPR = 1 - roc_rf$specificities,
  TPR = roc_rf$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#0072B2", size = 1.2) +
  geom_abline(linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  labs(
    title = "Random Forest ROC — TNBC vs non-TNBC",
    subtitle = paste0("AUC = ", round(auc_rf, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

ggsave(file.path(out_dir, "ROC_RF.pdf"), p_roc, width = 7, height = 6)
ggsave(file.path(out_dir, "ROC_RF.png"), p_roc, width = 7, height = 6, dpi = 300)

# -----------------------------
# Integrate everything for FINAL ranking
# -----------------------------
cat("🧠 Integrating ML + Network + DE + Mutation evidence...\n")

hub2 <- hub[, .(Gene, Degree)]
de2  <- de[, .(Gene, logFC, adj.P.Val)]
mut2 <- mut[, .(Gene)]

final <- data.frame(Gene = feature_genes)

final <- merge(final, rf_imp_df, by = "Gene", all.x = TRUE)
final <- merge(final, coef_df, by = "Gene", all.x = TRUE)
final <- merge(final, hub2, by = "Gene", all.x = TRUE)
final <- merge(final, de2, by = "Gene", all.x = TRUE)

final$Is_Mutated <- final$Gene %in% mut2$Gene

# Replace NAs
final[is.na(final)] <- 0

# Composite score
final$Score <- 
  scale(final$RF_Importance) +
  scale(abs(final$LASSO_Coefficient)) +
  scale(final$Degree) +
  scale(abs(final$logFC)) +
  ifelse(final$Is_Mutated, 1, 0)

final <- final[order(-final$Score), ]

# Save
fwrite(final, file.path(out_dir, "Final_target_ranking.tsv"), sep = "\t")

# -----------------------------
# Save summary
# -----------------------------
summary_file <- file.path(out_dir, "ML_QC_summary.txt")
sink(summary_file)

cat("STEP 9: Network-guided ML summary\n")
cat("=================================\n")
cat("Samples:", nrow(X), "\n")
cat("Features:", ncol(X), "\n")
cat("RF AUC:", auc_rf, "\n")
cat("Non-zero LASSO features:", nrow(coef_df), "\n")
cat("Final ranked genes:", nrow(final), "\n")

sink()

cat("✅ STEP 9 COMPLETE.\n")
cat("📁 Results saved to:", out_dir, "\n")
