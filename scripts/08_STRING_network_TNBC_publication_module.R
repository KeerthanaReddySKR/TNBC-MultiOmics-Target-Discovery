############################################
# STEP 8 — STRING Network & Hub Discovery
# Publication-grade Systems Biology Module
############################################

cat("🧬 Starting STEP 8: STRING network & hub discovery...\n")

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(STRINGdb)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
base_dir <- "D:/TNBC_Project_CLEAN"

de_file   <- file.path(base_dir, "results/DE_TNBC_vs_nonTNBC/DE_full_TNBC_vs_nonTNBC.tsv")
gsea_file <- file.path(base_dir, "results/GSEA_Hallmark/GSEA_Hallmark_full_results.tsv")
mut_file  <- file.path(base_dir, "results/MUT_TNBC_vs_nonTNBC/Mutation_comparison_significant.tsv")

out_dir <- file.path(base_dir, "results/STRING_Network")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Check inputs
# -----------------------------
stopifnot(file.exists(de_file))
stopifnot(file.exists(gsea_file))
stopifnot(file.exists(mut_file))

cat("✅ All input files found.\n")

# -----------------------------
# Load data
# -----------------------------
cat("📥 Loading DE results...\n")
de <- fread(de_file)

cat("📥 Loading GSEA results...\n")
gsea <- fread(gsea_file)

cat("📥 Loading mutation results...\n")
mut <- fread(mut_file)

# -----------------------------
# Sanity checks
# -----------------------------
if (!"Gene" %in% colnames(de)) {
  stop("DE file must contain a column named 'Gene'")
}

if (!"adj.P.Val" %in% colnames(de)) {
  stop("DE file must contain a column named 'adj.P.Val'")
}

if (!"logFC" %in% colnames(de)) {
  stop("DE file must contain a column named 'logFC'")
}

if (!"padj" %in% colnames(gsea)) {
  stop("GSEA file must contain a column named 'padj' (fgsea FDR column)")
}

# -----------------------------
# Select significant DE genes
# -----------------------------
de_sig <- de[ adj.P.Val < 0.05 & abs(logFC) >= 1 ]

cat("Significant DE genes:", nrow(de_sig), "\n")

# -----------------------------
# Select significant GSEA pathways (for QC)
# -----------------------------
gsea_sig <- gsea[ padj < 0.05 ]
cat("Significant Hallmark pathways:", nrow(gsea_sig), "\n")

# -----------------------------
# Mutation genes
# -----------------------------
mut_genes <- unique(mut$Gene)
cat("Significant mutation genes:", length(mut_genes), "\n")

# -----------------------------
# Build TNBC gene universe
# -----------------------------
tnbc_genes <- unique(c(
  de_sig$Gene,
  mut_genes
))

cat("Total TNBC gene universe size:", length(tnbc_genes), "\n")

write.table(
  data.frame(Gene = tnbc_genes),
  file = file.path(out_dir, "TNBC_gene_universe.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# Initialize STRING
# -----------------------------
cat("🌐 Initializing STRINGdb...\n")

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

# -----------------------------
# Map genes to STRING IDs
# -----------------------------
cat("🔗 Mapping genes to STRING IDs...\n")

map_df <- data.frame(Gene = tnbc_genes)
map_df <- string_db$map(map_df, "Gene", removeUnmappedRows = TRUE)

cat("Mapped genes to STRING:", nrow(map_df), "\n")

# -----------------------------
# Retrieve network
# -----------------------------
cat("🕸️ Retrieving STRING interactions...\n")

hits <- unique(map_df$STRING_id)
edges <- string_db$get_interactions(hits)

# Build igraph object
g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE)

cat("Network nodes (raw):", vcount(g), "\n")
cat("Network edges (raw):", ecount(g), "\n")

# Remove isolated nodes
g <- delete.vertices(g, degree(g) == 0)

cat("After pruning - nodes:", vcount(g), "\n")
cat("After pruning - edges:", ecount(g), "\n")

# -----------------------------
# Compute network metrics
# -----------------------------
cat("📊 Computing network centrality metrics...\n")

deg <- degree(g)
bet <- betweenness(g)
clo <- closeness(g)

net_table <- data.frame(
  STRING_id = names(deg),
  Degree = as.numeric(deg),
  Betweenness = as.numeric(bet),
  Closeness = as.numeric(clo)
)

# Add gene symbols
net_table <- merge(
  net_table,
  map_df[, c("STRING_id", "Gene")],
  by = "STRING_id"
)

# -----------------------------
# Rank hubs
# -----------------------------
net_table <- net_table[order(-net_table$Degree), ]

write.table(
  net_table,
  file = file.path(out_dir, "STRING_network_hub_table.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# Top 20 hubs
# -----------------------------
top_hubs <- head(net_table, 20)

write.table(
  top_hubs,
  file = file.path(out_dir, "Top20_STRING_hubs.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# Plot top hubs
# -----------------------------
p <- ggplot(top_hubs, aes(x = reorder(Gene, Degree), y = Degree)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 20 Hub Genes in TNBC Network (STRING)",
    x = "Gene",
    y = "Degree (Number of Interactions)"
  )

ggsave(
  filename = file.path(out_dir, "Top20_STRING_hubs.png"),
  plot = p,
  width = 7, height = 6, dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Top20_STRING_hubs.pdf"),
  plot = p,
  width = 7, height = 6
)

# -----------------------------
# Save network object
# -----------------------------
saveRDS(g, file = file.path(out_dir, "TNBC_STRING_network.rds"))

# -----------------------------
# QC report
# -----------------------------
qc_file <- file.path(out_dir, "STRING_QC_summary.txt")

qc_text <- c(
  "STEP 8 — STRING Network QC Summary",
  "----------------------------------",
  paste("DE significant genes:", nrow(de_sig)),
  paste("Mutation significant genes:", length(mut_genes)),
  paste("TNBC gene universe:", length(tnbc_genes)),
  paste("Mapped to STRING:", nrow(map_df)),
  paste("Network nodes:", vcount(g)),
  paste("Network edges:", ecount(g))
)

writeLines(qc_text, qc_file)

cat("✅ STEP 8 COMPLETE.\n")
cat("📁 Results saved to:", out_dir, "\n")
