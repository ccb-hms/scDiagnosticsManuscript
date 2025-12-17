# ----------------------------------------------
# COVID-19 PBMC - UMAP Visualization (FastMNN)
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
library(batchelor)      
library(BiocParallel)   
library(ggplot2)
library(dplyr)
library(pals) 
library(viridis)
library(scattermore) 

# ----------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("COVID-19 UMAP VISUALIZATION (FastMNN)\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading processed datasets...\n")
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")
cat(sprintf("Normal dataset: %d genes × %d cells\n", nrow(normal_data), ncol(normal_data)))
cat(sprintf("COVID dataset: %d genes × %d cells\n", nrow(covid_data), ncol(covid_data)))

# Setting the seed
set.seed(0)

# __________________________________________________
# Combine Datasets & Create Hybrid Cell Type Column
# __________________________________________________

cat("\nCombining datasets and creating hybrid cell type column...\n")

# Find common genes
common_genes <- intersect(rownames(normal_data), rownames(covid_data))
cat(sprintf("Found %d common genes.\n", length(common_genes)))

# Combine assays
combined_assay <- cbind(
    assay(normal_data, "logcounts")[common_genes, ],
    assay(covid_data, "logcounts")[common_genes, ]
)

# Create a unified colData with all necessary columns 
# For normal data, we use its own ground truth
normal_coldata_df <- as.data.frame(colData(normal_data)[, c("author_cell_type_merged", "Site")])
normal_coldata_df$disease_status <- "Healthy"
normal_coldata_df$hybrid_cell_type <- normal_coldata_df$author_cell_type_merged

# For COVID data, we include its Azimuth prediction
covid_coldata_df <- as.data.frame(colData(covid_data)[, c("author_cell_type_merged", "Site", "azimuth_celltype_l1_merged")])
covid_coldata_df$disease_status <- "COVID-19"
covid_coldata_df$hybrid_cell_type <- covid_coldata_df$azimuth_celltype_l1_merged

# Combine the data frames
combined_coldata <- bind_rows(normal_coldata_df, covid_coldata_df)

# Create the final combined SingleCellExperiment object
combined_sce <- SingleCellExperiment(
    assays = list(logcounts = combined_assay),
    colData = DataFrame(combined_coldata)
)
cat(sprintf("✓ Combined SCE object created with %d cells and a 'hybrid_cell_type' column.\n", ncol(combined_sce)))

# _____________________________
# FastMNN and UMAP Computation
# _____________________________

cat("\nComputing batch correction with FastMNN and UMAP...\n")
set.seed(123)

diagnostic_hvgs <- metadata(normal_data)$diagnostic_hvgs
available_hvgs <- intersect(diagnostic_hvgs, rownames(combined_sce))
cat(sprintf("Using %d/%d diagnostic HVGs for integration.\n", length(available_hvgs), length(diagnostic_hvgs)))

# Run fastMNN
cat("Applying FastMNN batch correction on 'Site'...\n")
mnn_output <- fastMNN(
    combined_sce,
    batch = combined_sce$Site,
    subset.row = available_hvgs,
    d = 50,
    BPPARAM = SerialParam()
)
cat("✓ FastMNN run completed.\n")

# Add corrected dimensions back to our main object
reducedDim(combined_sce, "MNN") <- reducedDim(mnn_output, "corrected")
cat("✓ Corrected dimensions ('MNN') added to the combined SCE object.\n")

# Compute UMAP
cat("Computing UMAP on batch-corrected data...\n")
combined_sce <- runUMAP(
    combined_sce,
    dimred = "MNN",
    n_dimred = 50,
    min_dist = 0.4
)
cat("✓ UMAP computed.\n")

# __________________________
# IFN Signature Computation
# __________________________

cat("\nComputing IFN signature scores...\n")
ifn_genes <- c("BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", 
               "IFI44L", "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", 
               "OAS1", "OAS2", "PARP9", "PLSCR1", "SAMD9", "SAMD9L", 
               "SP110", "STAT1", "TRIM22", "UBE2L6", "XAF1", "IRF7")
available_ifn_genes <- intersect(ifn_genes, rownames(combined_sce))

if(length(available_ifn_genes) > 0) {
    ifn_expression <- assay(combined_sce, "logcounts")[available_ifn_genes, , drop = FALSE]
    combined_sce$IFN_signature <- colMeans(ifn_expression)
    cat(sprintf("✓ IFN signature score computed using %d available genes.\n", length(available_ifn_genes)))
} else {
    combined_sce$IFN_signature <- 0
}

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

cell_type_colors <- c(
    
    # T Cell & ILC Lineage
    "CD4 T"       = "#8A2BE2",  
    "CD8 T"       = "#FF8C00",  
    "Treg"        = "#DA70D6",  
    "NKT"         = "#9932CC",  
    "MAIT"        = "#FF4500",  
    "gdT"         = "#FFA500",  
    "ILC"         = "#FF69B4",  
    
    # B Cell & Plasma Lineage
    "B cell"      = "#20B2AA",  
    "Plasmablast" = "#32CD32",  
    "Plasma cell" = "#006400",  

    # Myeloid & Monocyte Lineage
    "CD14 mono"   = "#008080",  
    "CD16 mono"   = "#1E90FF",  
    "Mono_prolif" = "#ADD8E6",  
    "DC1"         = "#A0522D",  
    "DC2"         = "#D2691E",  
    "DC3"         = "#F4A460",  
    "ASDC"        = "#8B4513",  
    "pDC"         = "#CD853F",  
    "DC_prolif"   = "#FFE4B5",  
    
    # NK Cell Lineage
    "NK_16hi"     = "#B22222",  
    "NK_56hi"     = "#8B4513",  
    "NK_prolif"   = "#FFC0CB",  
    
    # Other Types
    "HSPC"        = "#F0E68C",  
    "Platelets"   = "#BDB76B",  
    "RBC"         = "#E9967A"   
)

# Define the disease status palette
disease_colors <- c("Healthy" = "#85C1E9", "COVID-19" = "#F5B7B1")
cat("✓ Custom color palettes created.\n")

cat("✓ Custom color palettes created.\n")

# Create a reusable ggplot theme
publication_theme <- theme_classic(base_size = 12) +
    theme(
        plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.line = element_line(colour = "black", linewidth = 0.5)
    )

# Define point aesthetics
point_size <- 0.1
point_alpha <- 0.3

# Plot 1: UMAP by disease status
cat("Creating UMAP plot by disease status...\n")  
p_disease <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = disease_status)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = disease_colors) + 
    labs(color = "Condition", x = "UMAP 1", y = "UMAP 2") +
    publication_theme +
    guides(color = guide_legend(override.aes = list(size=5)))

# Plot 2: UMAP by hybrid cell types (with no legend)
cat("Creating UMAP plot by hybrid cell type...\n")
p_celltype <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = hybrid_cell_type), size = point_size, alpha = point_alpha + 0.1) +
    scale_color_manual(values = cell_type_colors, na.value = "grey80") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme  +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) 

# Plot 3: UMAP by IFN signature
cat("Creating UMAP plot by IFN signature...\n")
p_ifn <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = IFN_signature)) +
    geom_point(size = point_size, alpha = 0.5) +
    scale_color_viridis_c(name = "IFN Score", option = "plasma") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme

cat("✓ UMAP plots created.\n")

# _______________________
# Display and Save Plots
# _______________________

cat("\nDisplaying and saving plots in order: Disease, Cell Type, IFN Signature...\n")

# Display the plots in the R session
print(p_disease)
print(p_celltype)
print(p_ifn)

# Create the output directory if it doesn't already exist
dir.create("figures/covid", showWarnings = FALSE, recursive = TRUE)

# Save the plots with descriptive names
ggsave("figures/covid/fastmnn_umap_disease.png", p_disease, width = 10, height = 8, dpi = 300)
ggsave("figures/covid/fastmnn_umap_celltype.png", p_celltype, width = 10, height = 8, dpi = 300)
ggsave("figures/covid/fastmnn_umap_ifn.png", p_ifn, width = 10, height = 8, dpi = 300)

cat("✓ Plots saved to figures/covid/ directory with 'fastmnn_' prefix.\n")
cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("FastMNN UMAP VISUALIZATION COMPLETED\n")
cat(paste(rep("=", 45), collapse = ""), "\n")