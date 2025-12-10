# --------------------------------------------------------------------
# Differential Expression Analysis to Define an IAF Signature
# --------------------------------------------------------------------

# Load necessary libraries
library(SingleCellExperiment)
library(scran)
library(dplyr)

# --------------------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\n--- Loading Day 0 (Healthy) and Day 9 (DSS) Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
cat("✓ Datasets loaded.\n")

# Setting the seed
set.seed(0)

# _______________________________________
# Prepare Data for Pairwise DE Analysis
# _______________________________________

cat("\n--- Preparing data for specific pairwise DE analysis ---\n")

# Define the exact cell populations for comparison
healthy_stromal_types <- c("Fibro 1", "Fibro 2", "Fibro 5", "Fibro 6", "Fibro 7", "Fibro 12", "Fibro 13")
inflamed_stromal_types <- c("IAF 1", "IAF 2", "IAF 3", "IAF 4", "IAF 5", "Fibro 4")

# Subset the original objects to get the cells we need
healthy_cells_sce <- healthy_data[, healthy_data$tier2 %in% healthy_stromal_types]
inflamed_cells_sce <- dss9_data[, dss9_data$tier2 %in% inflamed_stromal_types]

cat(sprintf("Found %d healthy stromal cells in Day 0 data.\n", ncol(healthy_cells_sce)))
cat(sprintf("Found %d inflamed stromal cells in Day 9 data.\n", ncol(inflamed_cells_sce)))

# Extract only the necessary components: logcounts and cell names
common_genes <- intersect(rownames(healthy_cells_sce), rownames(inflamed_cells_sce))
healthy_logcounts <- assay(healthy_cells_sce, "logcounts")[common_genes, ]
inflamed_logcounts <- assay(inflamed_cells_sce, "logcounts")[common_genes, ]

# Manually construct a NEW, combined SingleCellExperiment object
combined_logcounts <- cbind(healthy_logcounts, inflamed_logcounts)

combined_coldata <- DataFrame(
  de_group = c(
    rep("Healthy_Stromal_D0", ncol(healthy_logcounts)),
    rep("Inflamed_Stromal_D9", ncol(inflamed_logcounts))
  ),
  row.names = c(colnames(healthy_logcounts), colnames(inflamed_logcounts))
)

# Create the final object for DE analysis
de_sce <- SingleCellExperiment(
    assays = list(logcounts = combined_logcounts),
    colData = combined_coldata
)
cat("✓ Created a clean, combined SCE object specifically for DE analysis.\n")


# _____________________________________
# Run Differential Expression Analysis
# _____________________________________

cat("\n--- Running scran::findMarkers to compare Day 9 IAFs vs. Day 0 Stromal Cells ---\n")

# Now we run findMarkers on our clean, purpose-built object
de_results <- findMarkers(
    de_sce,
    groups = de_sce$de_group,
    test.type = "wilcox",
    direction = "down",
    lfc = 0.5, 
    pval.type = "any"
)

inflamed_markers <- de_results[["Inflamed_Stromal_D9"]]
significant_inflamed_markers <- inflamed_markers[inflamed_markers$FDR < 0.05, ]

cat(sprintf("✓ Found %d significant upregulated genes in Inflamed Stromal cells (FDR < 0.05, LFC > 0.5).\n", 
            nrow(significant_inflamed_markers)))


# __________________________________
# Extract and Display the Top Genes
# __________________________________

cat("\n--- Top 25 Genes for Data-Driven 'Inflamed Fibroblast Signature' ---\n")

if (nrow(significant_inflamed_markers) > 0) {
    top_genes_df <- as.data.frame(significant_inflamed_markers) %>%
        arrange(desc(summary.AUC)) %>% # Rank by AUC effect size
        head(25)
    
    top_gene_list <- rownames(top_genes_df)

    cat("\n\n")
    cat('fibro_inflammation_genes <- c(\n    "', paste(top_gene_list, collapse = '",\n    "'), '"\n)\n', sep = "")

    cat("\n\n--- Detailed stats for top genes ---\n")
    print(top_genes_df)
} else {
    cat("No significant marker genes found with the current thresholds.\n")
}

