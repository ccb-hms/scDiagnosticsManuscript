# ----------------------------------------------------
# MERFISH Mouse Colon IBD - UMAP Visualization
# ----------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
library(Matrix) # Efficient handling of sparse matrices

# ----------------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("MERFISH UMAP VISUALIZATION (scVI Results)\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading processed and annotated datasets...\n")
normal_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data   <- readRDS("data/merfish/dss9_data.rds")

# Helper function to find UMAP coordinates
get_umap_coords <- function(sce_obj, obj_name) {
  possible_names <- c("UMAP_scVI", "X_umap", "UMAP")
  found_name <- intersect(possible_names, reducedDimNames(sce_obj))[1]
  
  if (is.na(found_name)) {
    stop(paste("Error: No UMAP coordinates found in", obj_name))
  }
  
  cat(paste("✓ Found UMAP coordinates in", obj_name, "under slot:", found_name, "\n"))
  coords <- as.data.frame(reducedDim(sce_obj, found_name))
  colnames(coords)[1:2] <- c("UMAP1", "UMAP2")
  return(coords)
}

# Extract coordinates
umap_normal <- get_umap_coords(normal_data, "Healthy Data")
umap_dss9   <- get_umap_coords(dss9_data, "DSS9 Data")

# _____________________________________________
# Create Combined Plotting Data Frame
# _____________________________________________

cat("\nCombining data for plotting...\n")

# Check for required column
if (!"tier2_merged" %in% colnames(colData(normal_data)) || !"tier2_merged" %in% colnames(colData(dss9_data))) {
  stop("Column 'tier2_merged' missing. Run the merging script first.")
}

# Extract metadata
meta_normal <- as.data.frame(colData(normal_data)) |>
    select(tier2_merged) |>
    mutate(disease_status = "Healthy",
           cell_type = tier2_merged,
           barcode = rownames(colData(normal_data))) # Keep barcode for matching

meta_dss9 <- as.data.frame(colData(dss9_data)) |>
    select(tier2_merged) |>
    mutate(disease_status = "DSS9 (Colitis)",
           cell_type = tier2_merged,
           barcode = rownames(colData(dss9_data)))

# Combine DF
plot_df <- bind_rows(
    bind_cols(umap_normal, meta_normal),
    bind_cols(umap_dss9, meta_dss9)
)
cat(sprintf("✓ Combined plotting data frame created with %d cells.\n", nrow(plot_df)))

# __________________________
# ECM Score Computation
# __________________________

cat("\nComputing ECM Homeostasis Score...\n")

# 1. Define Signature
ecm_genes <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")

# 2. Check Gene Availability (MERFISH panels are small)
common_genes <- intersect(rownames(normal_data), rownames(dss9_data))
available_ecm_genes <- intersect(ecm_genes, common_genes)

cat("Requested genes:", paste(ecm_genes, collapse=", "), "\n")
cat("Available genes:", paste(available_ecm_genes, collapse=", "), "\n")

if(length(available_ecm_genes) == 0) {
    warning("None of the ECM genes were found in the dataset! Setting score to 0.")
    plot_df$ECM_Score <- 0
} else {
    # 3. Extract Expression Data
    # We extract data only for the available genes
    # Note: Using 'logcounts' as per your reference script
    
    expr_normal <- assay(normal_data, "logcounts")[available_ecm_genes, , drop=FALSE]
    expr_dss9   <- assay(dss9_data, "logcounts")[available_ecm_genes, , drop=FALSE]
    
    # 4. Calculate Raw Mean per cell
    scores_normal <- colMeans(as.matrix(expr_normal), na.rm = TRUE)
    scores_dss9   <- colMeans(as.matrix(expr_dss9), na.rm = TRUE)
    
    # Combine scores into a named vector
    all_scores <- c(scores_normal, scores_dss9)
    
    # 5. Apply Logic: Cap at 1.5 -> Multiply by 2
    # Cap (Winsorize)
    capped_scores <- pmin(all_scores, 1.5)
    
    # Scale
    final_scores <- capped_scores * 2
    
    # 6. Map back to plot_df
    # We use match() on barcodes to ensure perfect alignment
    plot_df$ECM_Score <- final_scores[match(plot_df$barcode, names(final_scores))]
    
    cat("✓ ECM scores computed (Range: ", min(plot_df$ECM_Score), " - ", max(plot_df$ECM_Score), ")\n")
}

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

# Define Color Palette (Provided)
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

# Disease Palette
disease_colors <- c("Healthy" = "#5DADE2", "DSS9 (Colitis)" = "#FF6347") 

# Theme
publication_theme <- theme_classic(base_size = 12) +
    theme(
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
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
    guides(color = guide_legend(override.aes = list(size = 5)))

# Plot 2: Cell Type
cat("Creating UMAP plot by cell type...\n")
types_to_keep <- plot_df |> count(cell_type) |> filter(n >= 10) |> pull(cell_type)
plot_df_filtered <- plot_df |> filter(cell_type %in% types_to_keep)

p_celltype <- ggplot(plot_df_filtered, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = cell_type), size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) 

# Plot 3: ECM Score
cat("Creating UMAP plot by ECM Score...\n")
p_ecm <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = ECM_Score)) +
    # Use slightly higher alpha for feature plots to see the gradient better
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

ggsave("figures/merfish/scvi_umap__disease.png", p_disease, width = 8, height = 6, dpi = 600)
ggsave("figures/merfish/scvi_umap__celltype.png", p_celltype, width = 10, height = 8, dpi = 600)
ggsave("figures/merfish/scvi_umap__ecm_score.png", p_ecm, width = 10, height = 8, dpi = 600)

cat("✓ Plots saved to figures/merfish/ directory.\n")
cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("scVI UMAP VISUALIZATION COMPLETED\n")
cat(paste(rep("=", 45), collapse = ""), "\n")