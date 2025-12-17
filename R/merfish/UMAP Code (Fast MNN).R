# ----------------------------------------------
# MERFISH Mouse Colon IBD - UMAP (FastMNN)
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
library(batchelor)      
library(BiocParallel)   
library(ggplot2)
library(dplyr)
library(viridis)
library(Matrix)

# ----------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("MERFISH UMAP VISUALIZATION (FastMNN)\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading processed datasets...\n")
normal_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data   <- readRDS("data/merfish/dss9_data.rds")

cat(sprintf("Healthy dataset: %d genes × %d cells\n", nrow(normal_data), ncol(normal_data)))
cat(sprintf("DSS9 dataset:    %d genes × %d cells\n", nrow(dss9_data), ncol(dss9_data)))

# Setting the seed
set.seed(0)

# __________________________________________________
# Combine Datasets & Create Hybrid Cell Type Column
# __________________________________________________

cat("\nCombining datasets...\n")

# Check for required columns
if (!"tier2_merged" %in% colnames(colData(normal_data))) stop("tier2_merged missing in Healthy data")
if (!"tier2_merged" %in% colnames(colData(dss9_data))) stop("tier2_merged missing in DSS9 data")

# Find common genes
common_genes <- intersect(rownames(normal_data), rownames(dss9_data))
cat(sprintf("Found %d common genes.\n", length(common_genes)))

# Combine assays (Logcounts)
combined_assay <- cbind(
    assay(normal_data, "logcounts")[common_genes, ],
    assay(dss9_data, "logcounts")[common_genes, ]
)

# Prepare Metadata
# Reference (Healthy)
normal_coldata_df <- as.data.frame(colData(normal_data)[, "tier2_merged", drop=FALSE])
normal_coldata_df$disease_status <- "Healthy"
normal_coldata_df$cell_type <- normal_coldata_df$tier2_merged
colnames(normal_coldata_df)[1] <- "original_annotation"

# Query (DSS9)
dss9_coldata_df <- as.data.frame(colData(dss9_data)[, "tier2_merged", drop=FALSE])
dss9_coldata_df$disease_status <- "DSS9 (Colitis)"
dss9_coldata_df$cell_type <- dss9_coldata_df$tier2_merged
colnames(dss9_coldata_df)[1] <- "original_annotation"

# Combine Metadata
combined_coldata <- bind_rows(normal_coldata_df, dss9_coldata_df)

# Restore rownames to match the combined_assay columns
rownames(combined_coldata) <- combined_coldata$barcode

# Create Combined SCE
combined_sce <- SingleCellExperiment(
    assays = list(logcounts = combined_assay),
    colData = DataFrame(combined_coldata)
)
cat(sprintf("✓ Combined SCE object created with %d cells.\n", ncol(combined_sce)))

# _____________________________
# FastMNN and UMAP Computation
# _____________________________

cat("\nComputing batch correction with FastMNN and UMAP...\n")
set.seed(123)

# For MERFISH, the gene panel is already selected for importance.
# We generally use all common genes for integration unless the panel is huge.
integration_genes <- rownames(combined_sce)
cat(sprintf("Using all %d common genes for MNN integration.\n", length(integration_genes)))

# Run fastMNN
# We treat 'disease_status' as the batch to align the datasets
cat("Applying FastMNN batch correction on 'disease_status'...\n")
mnn_output <- fastMNN(
    combined_sce,
    batch = combined_sce$disease_status,
    subset.row = integration_genes,
    d = 50,
    BPPARAM = SerialParam()
)
cat("✓ FastMNN run completed.\n")

# Add corrected dimensions back to our main object
reducedDim(combined_sce, "MNN") <- reducedDim(mnn_output, "corrected")

# Compute UMAP on MNN embeddings
cat("Computing UMAP on batch-corrected data...\n")
combined_sce <- runUMAP(
    combined_sce,
    dimred = "MNN",
    n_dimred = 50,
    min_dist = 0.4
)
cat("✓ UMAP computed.\n")

# __________________________
# ECM Score Computation
# __________________________

cat("\nComputing ECM Homeostasis Score...\n")

# 1. Define Signature
ecm_genes <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")
available_ecm_genes <- intersect(ecm_genes, rownames(combined_sce))

if(length(available_ecm_genes) > 0) {
    # 2. Extract Expression
    expr_matrix <- assay(combined_sce, "logcounts")[available_ecm_genes, , drop = FALSE]
    
    # 3. Raw Mean
    raw_score <- colMeans(as.matrix(expr_matrix), na.rm = TRUE)
    
    # 4. Cap at 1.5
    capped_score <- pmin(raw_score, 1.5)
    
    # 5. Scale by 2
    final_score <- capped_score * 2
    
    combined_sce$ECM_Score <- final_score
    cat(sprintf("✓ ECM score computed using %d genes.\n", length(available_ecm_genes)))
} else {
    combined_sce$ECM_Score <- 0
    warning("No ECM genes found in dataset.")
}

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

# Extract plotting data
plot_df <- as.data.frame(reducedDim(combined_sce, "UMAP"))
colnames(plot_df) <- c("UMAP1", "UMAP2")
plot_df <- cbind(plot_df, as.data.frame(colData(combined_sce)))

# Define Colors
bg_color <- "grey90"
merfish_colors <- c(
    "Inflamed Fibroblast" = "#A50F15",  
    "Fibroblast"          = "#FB9A99",  
    "Inflamed SMC"        = "#54278F",  
    "Smooth Muscle"       = "#BC80BD",  
    "Inflamed Epithelial" = "#E7298A",  
    "Epithelial"          = "#F1B6DA",  
    "Stem/TA"             = "#008B8B",  
    "Other Immune"        = "#FDB462",  
    "Neutrophil"          = "#4682B4",  
    "Pericyte/FRC"        = "#8C510A",  
    "Endothelial"         = "#B2DF8A",  
    "Enteric Nervous"     = "#FFFF99",  
    "ICC"                 = "#8DD3C7",  
    "Adipose"             = "#E5C494",  
    "Other"               = bg_color
)

disease_colors <- c("Healthy" = "#5DADE2", "DSS9 (Colitis)" = "#FF6347") 

# Theme
publication_theme <- theme_classic(base_size = 12) +
    theme(
        plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.line = element_line(colour = "black", linewidth = 0.5)
    )

point_size <- 0.1
point_alpha <- 0.6

# Plot 1: Disease Status
cat("Creating UMAP plot by disease status...\n")  
p_disease <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = disease_status)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = disease_colors) + 
    labs(color = "Condition", x = "UMAP 1", y = "UMAP 2") +
    publication_theme +
    guides(color = guide_legend(override.aes = list(size=5)))

# Plot 2: Cell Type
cat("Creating UMAP plot by cell type...\n")

# Filter legends slightly
types_to_keep <- plot_df |> count(cell_type) |> filter(n >= 10) |> pull(cell_type)
plot_df_filtered <- plot_df |> filter(cell_type %in% types_to_keep)

p_celltype <- ggplot(plot_df_filtered, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = cell_type), size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme  +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) 

# Plot 3: ECM Score
cat("Creating UMAP plot by ECM Score...\n")
p_ecm <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = ECM_Score)) +
    geom_point(size = point_size, alpha = 0.8) +
    scale_color_viridis_c(name = "ECM Score", option = "magma") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme

cat("✓ UMAP plots created.\n")

# _______________________
# Display and Save Plots
# _______________________

cat("\nDisplaying and saving plots...\n")

print(p_disease)
print(p_celltype)
print(p_ecm)

dir.create("figures/merfish", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/merfish/fastmnn_umap_disease.png", p_disease, width = 8, height = 6, dpi = 600)
ggsave("figures/merfish/fastmnn_umap_celltype.png", p_celltype, width = 10, height = 8, dpi = 600)
ggsave("figures/merfish/fastmnn_umap_ecm.png", p_ecm, width = 10, height = 8, dpi = 600)

cat("✓ Plots saved to figures/merfish/ directory with 'fastmnn_' prefix.\n")
cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("FastMNN UMAP VISUALIZATION COMPLETED\n")
cat(paste(rep("=", 45), collapse = ""), "\n")