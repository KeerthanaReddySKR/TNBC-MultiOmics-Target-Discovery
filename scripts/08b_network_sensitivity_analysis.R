# ================================
# STEP 8B — Network Sensitivity / Robustness Analysis
# ================================

rm(list = ls())
cat("🧪 Starting STEP 8B: Network sensitivity analysis...\n")

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(STRINGdb)
})

base_dir <- "D:/TNBC_Project_CLEAN"

ranking_file <- file.path(base_dir, "results/ML_Target_Ranking/Final_target_ranking.tsv")
out_dir <- file.path(base_dir, "results/STRING_Network/Sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(ranking_file))

# Load ranking
rank <- fread(ranking_file)

# Define cutoffs
cutoffs <- c(20, 30, 40)

# Initialize STRING
cat("🌐 Initializing STRINGdb...\n")
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

hub_lists <- list()

for (k in cutoffs) {
  
  cat("🔬 Processing Top", k, "genes...\n")
  
  genes_k <- rank$Gene[1:k]
  
  # Map to STRING
  map_df <- data.frame(Gene = genes_k)
  mapped <- string_db$map(map_df, "Gene", removeUnmappedRows = TRUE)
  
  # Get interactions
  edges <- string_db$get_interactions(mapped$STRING_id)
  
  if (nrow(edges) == 0) {
    cat("⚠️ No edges found for Top", k, "\n")
    next
  }
  
  # Build graph
  g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE)
  
  # Compute centralities
  deg <- degree(g)
  bet <- betweenness(g)
  
  metrics <- data.frame(
    STRING_id = names(deg),
    Degree = as.numeric(deg),
    Betweenness = as.numeric(bet)
  )
  
  # Map back to gene symbols
  inv_map <- mapped[, c("STRING_id", "Gene")]
  metrics <- merge(metrics, inv_map, by = "STRING_id")
  
  metrics <- metrics[order(-metrics$Degree), ]
  
  # Save
  out_file <- file.path(out_dir, paste0("Top", k, "_network_metrics.tsv"))
  fwrite(metrics, out_file, sep = "\t")
  
  # Store top hubs
  top_hubs <- head(metrics$Gene, 10)
  hub_lists[[paste0("Top", k)]] <- top_hubs
  
  cat("✅ Top hubs for Top", k, ":", paste(top_hubs, collapse = ", "), "\n")
}

# Compare overlaps
cat("📊 Computing overlap...\n")

all_hubs <- unique(unlist(hub_lists))

overlap_df <- data.frame(
  Gene = all_hubs,
  In_Top20 = all_hubs %in% hub_lists$Top20,
  In_Top30 = all_hubs %in% hub_lists$Top30,
  In_Top40 = all_hubs %in% hub_lists$Top40
)

fwrite(overlap_df, file.path(out_dir, "Hub_overlap_summary.tsv"), sep = "\t")

cat("✅ Sensitivity analysis complete.\n")
cat("📁 Results in:", out_dir, "\n")
