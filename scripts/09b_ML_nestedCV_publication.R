############################################
# STEP 9B — Nested Cross-Validation ML
# Publication-grade, reviewer-safe
############################################

cat("🧠 Starting STEP 9B: Nested CV ML (Reviewer-grade)...\n")

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(caret)
  library(data.table)
  library(ggplot2)
})

set.seed(123)

# -----------------------------
# Paths
# -----------------------------
base_dir <- "D:/TNBC_Project_CLEAN"

expr_file <- file.path(base_dir, "data_processed/expr_1218.rds")
feature_file <- file.path(base_dir, "results/ML_Target_Ranking/Final_target_ranking.tsv")

out_dir <- file.path(base_dir, "results/ML_NestedCV")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------
cat("📥 Loading data...\n")

obj <- readRDS(expr_file)

expr <- obj$expr
labels <- obj$labels

features <- fread(feature_file)

# -----------------------------
# Use EXISTING TNBC label
# -----------------------------
stopifnot("is_TNBC" %in% colnames(labels))

y_raw <- labels$is_TNBC

# 🔴 CRITICAL FIX: convert logical -> 0/1
if (is.logical(y_raw)) {
  y_raw <- ifelse(y_raw, 1, 0)
}

# Safety check
if (!all(y_raw %in% c(0,1,NA))) {
  stop("❌ is_TNBC must be 0/1 or TRUE/FALSE")
}

cat("TNBC samples:", sum(y_raw == 1, na.rm = TRUE), "\n")
cat("Non-TNBC samples:", sum(y_raw == 0, na.rm = TRUE), "\n")

# -----------------------------
# Select top features
# -----------------------------
top_features <- head(features$Gene, 200)
top_features <- intersect(top_features, rownames(expr))

cat("Features used:", length(top_features), "\n")

# -----------------------------
# Build ML matrix
# -----------------------------
X <- t(expr[top_features, ])
X <- scale(X)

# -----------------------------
# Filter samples
# -----------------------------
keep <- !is.na(y_raw)

X <- X[keep, ]
y_raw <- y_raw[keep]

y <- factor(y_raw, levels = c(0,1), labels = c("nonTNBC","TNBC"))

cat("Samples after filtering:", nrow(X), "\n")
cat("Features:", ncol(X), "\n")
print(table(y))

# -----------------------------
# Safety check
# -----------------------------
if (any(table(y) < 10)) {
  stop("❌ One class has too few samples. Check labels.")
}

# -----------------------------
# Nested CV setup
# -----------------------------
outer_folds <- createFolds(y, k = 10, returnTrain = FALSE)

auc_outer <- c()
roc_list <- list()

fold_id <- 1

# -----------------------------
# Outer loop
# -----------------------------
for (test_idx in outer_folds) {
  
  cat("🔁 Outer fold:", fold_id, "\n")
  
  X_train <- X[-test_idx, ]
  y_train <- y[-test_idx]
  
  X_test <- X[test_idx, ]
  y_test <- y[test_idx]
  
  # -----------------------------
  # Inner CV (5-fold) for lambda
  # -----------------------------
  cvfit <- cv.glmnet(
    x = as.matrix(X_train),
    y = y_train,
    family = "binomial",
    alpha = 1,
    nfolds = 5,
    type.measure = "auc"
  )
  
  best_lambda <- cvfit$lambda.min
  
  # -----------------------------
  # Train final model
  # -----------------------------
  model <- glmnet(
    x = as.matrix(X_train),
    y = y_train,
    family = "binomial",
    alpha = 1,
    lambda = best_lambda
  )
  
  # -----------------------------
  # Predict on test
  # -----------------------------
  prob <- predict(model, newx = as.matrix(X_test), type = "response")[,1]
  
  roc_obj <- roc(y_test, prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_outer <- c(auc_outer, auc_val)
  roc_list[[fold_id]] <- roc_obj
  
  cat("  Fold AUC:", round(auc_val, 3), "\n")
  
  fold_id <- fold_id + 1
}

# -----------------------------
# Results
# -----------------------------
mean_auc <- mean(auc_outer)
sd_auc <- sd(auc_outer)

cat("=====================================\n")
cat("Nested CV AUC (mean ± sd):", round(mean_auc,3), "±", round(sd_auc,3), "\n")
cat("=====================================\n")

# -----------------------------
# Save QC summary
# -----------------------------
qc_text <- c(
  "STEP 9B — Nested Cross-Validation ML QC",
  "--------------------------------------",
  paste("Folds: 10 outer x 5 inner"),
  paste("Features used:", ncol(X)),
  paste("Samples used:", nrow(X)),
  paste("TNBC:", sum(y == "TNBC")),
  paste("Non-TNBC:", sum(y == "nonTNBC")),
  paste("Mean AUC:", round(mean_auc, 4)),
  paste("SD AUC:", round(sd_auc, 4)),
  "",
  "Per-fold AUC:",
  paste(round(auc_outer, 4), collapse = ", ")
)

writeLines(qc_text, file.path(out_dir, "ML_NestedCV_QC_summary.txt"))

# -----------------------------
# Plot ROC curves
# -----------------------------
png(file.path(out_dir, "NestedCV_ROC.png"), width = 2000, height = 1600, res = 200)

plot(roc_list[[1]], col = "blue", lwd = 2, main = "Nested CV ROC — TNBC vs non-TNBC")
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = rgb(0,0,1,0.3), lwd = 1, add = TRUE)
}
abline(a=0,b=1,lty=2,col="gray")

legend("bottomright", legend = paste0("Mean AUC = ", round(mean_auc,3)), bty = "n")

dev.off()

cat("✅ STEP 9B COMPLETE. Results in:", out_dir, "\n")
