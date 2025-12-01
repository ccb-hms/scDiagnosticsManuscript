# --------------------------------------------------------------
# COVID-19 - Dataset Combination and UMAP Visualization (scVI)
# --------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(reticulate)
library(zellkonverter)

# --------------------------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("COVID-19 DATASET COMBINATION AND UMAP VISUALIZATION\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading processed datasets...\n")

# Load both processed datasets
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

cat(sprintf("Normal dataset: %d genes × %d cells\n", nrow(normal_data), ncol(normal_data)))
cat(sprintf("COVID dataset: %d genes × %d cells\n", nrow(covid_data), ncol(covid_data)))

# Verify required columns exist
required_cols <- c("author_cell_type", 
                   "author_cell_type_merged", 
                   "Status_on_day_collection_summary",
                   "Site", 
                   "sample_id")
for(col in required_cols) {
    if(!col %in% colnames(colData(normal_data))) {
        stop(sprintf("Column '%s' not found in normal dataset", col))
    }
    if(!col %in% colnames(colData(covid_data))) {
        stop(sprintf("Column '%s' not found in COVID dataset", col))
    }
}

cat("✓ Required columns verified in both datasets\n")

# ___________________________________
# Combine Datasets for Visualization
# ___________________________________

cat("\nCombining normal and COVID datasets for visualization...\n")

# Extract assays (logcounts)
normal_assay <- assay(normal_data, "logcounts")
covid_assay <- assay(covid_data, "logcounts")

cat(sprintf("Normal assay dimensions: %d genes × %d cells\n", nrow(normal_assay), ncol(normal_assay)))
cat(sprintf("COVID assay dimensions: %d genes × %d cells\n", nrow(covid_assay), ncol(covid_assay)))

# Ensure same genes in same order
common_genes <- intersect(rownames(normal_assay), rownames(covid_assay))
normal_assay <- normal_assay[common_genes, ]
covid_assay <- covid_assay[common_genes, ]

cat(sprintf("Common genes: %d\n", length(common_genes)))

# Combine assays
combined_assay <- cbind(normal_assay, covid_assay)
cat(sprintf("Combined assay dimensions: %d genes × %d cells\n", nrow(combined_assay), ncol(combined_assay)))

# Extract relevant colData columns
normal_coldata <- colData(normal_data)[, c("author_cell_type", 
                                           "author_cell_type_merged", 
                                           "Status_on_day_collection_summary",
                                           "Site", 
                                           "sample_id")]
covid_coldata <- colData(covid_data)[, c("author_cell_type", 
                                         "author_cell_type_merged", 
                                         "Status_on_day_collection_summary",
                                         "Site", 
                                         "sample_id")]

# Add disease status column
normal_coldata$disease_status <- "Healthy"
covid_coldata$disease_status <- "COVID-19"

# Combine colData
combined_coldata <- rbind(normal_coldata, covid_coldata)

cat(sprintf("Combined colData dimensions: %d cells × %d columns\n", nrow(combined_coldata), ncol(combined_coldata)))

# Create combined SCE object
combined_sce <- SingleCellExperiment(
    assays = list(logcounts = combined_assay),
    colData = combined_coldata
)

cat("✓ Combined SCE object created\n")

# __________________________
# scVI and UMAP Computation
# __________________________

cat("\nSetting up Python environment via Reticulate...\n")
# Point reticulate to your conda environment
# Make sure the environment name 'scvi-env' matches what you created
tryCatch({
    use_condaenv("scvi-env", required = TRUE)
    cat("✓ Python environment 'scvi-env' loaded.\n")
}, error = function(e) {
    stop("Could not find or load the 'scvi-env' conda environment. Please ensure it's created and has scvi-tools installed.")
})

cat("\nImporting required Python libraries...\n")
# Import python libraries into your R session
sc <- import("scanpy")
scvi <- import("scvi")
np <- import("numpy")

cat("✓ Python libraries imported.\n")

cat("\nConverting R SingleCellExperiment to Python AnnData object...\n")

# This step can consume a lot of memory.
adata <- SCE2AnnData(combined_sce, main_layer = "counts")
cat("✓ Conversion to AnnData complete.\n")

cat("\nSetting up and training the scVI model...\n")
# This is the equivalent of the Python workflow
scvi$model$SCVI$setup_anndata(adata, batch_key = 'Site')

# Create the model
# n_latent=30 is a common choice, but can be tuned
vae <- scvi$model$SCVI(adata, n_latent = 30)

# Train the model. 
# It will benefit greatly from a GPU if scvi-tools was installed with GPU support.
# On a CPU, this can take hours for a large dataset.
vae$train()
cat("✓ scVI model training complete.\n")

cat("\nExtracting latent space and running UMAP...\n")
# Get the batch-corrected latent representation
latent_representation <- vae$get_latent_representation()

# Add it to the AnnData object's .obsm slot
adata$obsm[["X_scVI"]] <- latent_representation

# Build the k-NN graph with the paper's parameters
sc$pp$neighbors(adata, use_rep = 'X_scVI', n_neighbors = 100L) # Use 100L for integer

# Run UMAP on the k-NN graph
sc$tl$umap(adata)
cat("✓ UMAP computed on scVI latent space.\n")

cat("\nTransferring UMAP results back to R...\n")
# Extract the UMAP coordinates from the AnnData object
umap_coords_scvi <- adata$obsm[['X_umap']]

# Add the new UMAP coordinates to your original R SCE object
# We name it "UMAP_scVI" to distinguish it from other UMAPs.
reducedDim(combined_sce, "UMAP_scVI") <- umap_coords_scvi

cat("✓ scVI workflow complete. UMAP coordinates added to SCE object.\n")

# _____________________________
# IFN Signature Computation
# _____________________________

cat("\nComputing IFN signature scores...\n")

# Define IFN-associated genes from Yoshida et al.
ifn_genes <- c("BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", 
               "IFI44L", "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", 
               "OAS1", "OAS2", "PARP9", "PLSCR1", "SAMD9", "SAMD9L", 
               "SP110", "STAT1", "TRIM22", "UBE2L6", "XAF1", "IRF7")

cat(sprintf("Total IFN signature genes: %d\n", length(ifn_genes)))

# Check which IFN genes are present in the dataset
available_ifn_genes <- intersect(ifn_genes, rownames(combined_sce))
missing_ifn_genes <- setdiff(ifn_genes, rownames(combined_sce))

cat(sprintf("Available IFN genes: %d/%d\n", length(available_ifn_genes), length(ifn_genes)))
if(length(missing_ifn_genes) > 0) {
    cat("Missing IFN genes:", paste(missing_ifn_genes, collapse = ", "), "\n")
}

# Compute IFN signature score (mean expression of available IFN genes)
if(length(available_ifn_genes) > 0) {
    
    # Extract expression matrix for IFN genes
    ifn_expression <- assay(combined_sce, "logcounts")[available_ifn_genes, , drop = FALSE]
    
    # Calculate mean expression across IFN genes for each cell
    ifn_signature_scores <- colMeans(ifn_expression)
    
    # Add IFN signature to colData
    combined_sce$IFN_signature <- ifn_signature_scores
    
    cat("✓ IFN signature scores computed and added to colData\n")
    
    # Summary statistics
    cat(sprintf("IFN signature range: %.3f to %.3f\n", 
        min(ifn_signature_scores), max(ifn_signature_scores)))
    cat(sprintf("Mean IFN signature: %.3f\n", mean(ifn_signature_scores)))
    
    # Compare IFN signature between healthy and COVID
    healthy_ifn <- ifn_signature_scores[combined_sce$disease_status == "Healthy"]
    covid_ifn <- ifn_signature_scores[combined_sce$disease_status == "COVID-19"]
    
    cat(sprintf("Mean IFN signature - Healthy: %.3f\n", mean(healthy_ifn)))
    cat(sprintf("Mean IFN signature - COVID-19: %.3f\n", mean(covid_ifn)))
    cat(sprintf("Fold change (COVID/Healthy): %.2f\n", mean(covid_ifn) / mean(healthy_ifn)))
    
    # Store IFN gene information in metadata
    metadata(combined_sce)$ifn_genes_all <- ifn_genes
    metadata(combined_sce)$ifn_genes_available <- available_ifn_genes
    metadata(combined_sce)$ifn_genes_missing <- missing_ifn_genes
    
} else {
    cat("Warning: No IFN signature genes found in dataset\n")
    combined_sce$IFN_signature <- rep(0, ncol(combined_sce))
}

cat("✓ IFN signature computation completed\n")

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

# Set up color palettes
# Color palette for cell types
n_cell_types <- length(unique(combined_sce$author_cell_type_merged))
cell_type_colors <- RColorBrewer::brewer.pal(min(n_cell_types, 12), "Set3")
if(n_cell_types > 12) {
    # Extend palette if more than 12 cell types
    cell_type_colors <- c(cell_type_colors, rainbow(n_cell_types - 12))
}
names(cell_type_colors) <- sort(unique(combined_sce$author_cell_type_merged))

# Color palette for disease status
disease_colors <- c("Healthy" = "#2E86C1", "COVID-19" = "#E74C3C")

# Store colors in metadata for consistency
metadata(combined_sce)$colors_cell_type <- cell_type_colors
metadata(combined_sce)$colors_disease <- disease_colors

# Plot 1: UMAP colored by merged cell types
cat("Creating UMAP plot colored by cell type...\n")
p1 <- plotReducedDim(combined_sce, "UMAP", colour_by = "author_cell_type_merged", 
                     scattermore = TRUE, rasterise = TRUE) +
    scale_color_manual(values = metadata(combined_sce)$colors_cell_type) + 
    labs(color = "Cell Type", 
         title = "UMAP: COVID-19 vs Healthy PBMCs (Corrected with FastMNN)",
         subtitle = "Colored by Cell Type") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

# Plot 2: UMAP colored by disease status
cat("Creating UMAP plot colored by disease status...\n")  
p2 <- plotReducedDim(combined_sce, "UMAP", colour_by = "disease_status", 
                     scattermore = TRUE, rasterise = TRUE) +
    scale_color_manual(values = metadata(combined_sce)$colors_disease) + 
    labs(color = "Disease Status", 
         title = "UMAP: COVID-19 vs Healthy PBMCs (Corrected with FastMNN)",
         subtitle = "Colored by Disease Status") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

# Plot 3: UMAP colored by IFN signature
cat("Creating UMAP plot colored by IFN signature...\n")
p3 <- plotReducedDim(combined_sce, "UMAP", colour_by = "IFN_signature", 
                     scattermore = TRUE, rasterise = TRUE) +
    scale_color_viridis_c(name = "IFN Signature", option = "plasma") + 
    labs(color = "IFN Signature", 
         title = "UMAP: COVID-19 vs Healthy PBMCs (Corrected with FastMNN)",
         subtitle = "Colored by IFN Signature Score") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

cat("✓ UMAP plots created\n")

# _______________________
# Display and Save Plots
# _______________________

cat("\nDisplaying plots...\n")

# Display plots
print(p1)
print(p2)
print(p3)

# Save plots
cat("Saving plots...\n")
ggsave("figures/covid/covid_umap_by_celltype_fastmnn.png", p1, width = 12, height = 8, dpi = 600)
ggsave("figures/covid/covid_umap_by_disease_fastmnn.png", p2, width = 12, height = 8, dpi = 600)
ggsave("figures/covid/covid_umap_by_ifn_signature_fastmnn.png", p3, width = 12, height = 8, dpi = 600)

cat("✓ Plots saved to figures/covid/ directory\n")

# Print summary statistics
cat("\nDataset summary:\n")
cat(sprintf("Total cells: %d\n", ncol(combined_sce)))
cat(sprintf("Healthy cells: %d (%.1f%%)\n", 
    sum(combined_sce$disease_status == "Healthy"),
    100 * sum(combined_sce$disease_status == "Healthy") / ncol(combined_sce)))
cat(sprintf("COVID-19 cells: %d (%.1f%%)\n", 
    sum(combined_sce$disease_status == "COVID-19"),
    100 * sum(combined_sce$disease_status == "COVID-19") / ncol(combined_sce)))
cat(sprintf("Number of cell types: %d\n", length(unique(combined_sce$author_cell_types_merged))))
cat(sprintf("IFN signature range: %.3f - %.3f\n", 
    min(combined_sce$IFN_signature), max(combined_sce$IFN_signature)))

cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("UMAP VISUALIZATION COMPLETED\n")
cat(paste(rep("=", 45), collapse = ""), "\n")