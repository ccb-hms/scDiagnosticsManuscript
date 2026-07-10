# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 2 Code
# -------------------------------------------------------

library(SingleCellExperiment)
library(SpatialExperiment) 
library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)
library(Matrix)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\nLoading DSS9 dataset for annotation validation...\n")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Robust helper function to find spatial coordinates anywhere in the object
get_spatial_coords <- function(sce_obj, obj_name) {
  
  # Method 1: Try SpatialExperiment::spatialCoords
  tryCatch({
    coords <- as.data.frame(spatialCoords(sce_obj))
    if (nrow(coords) > 0) {
      cat(paste("✓ Found spatial coordinates in", obj_name, "via spatialCoords()\n"))
      colnames(coords)[1:2] <- c("x", "y")
      return(coords)
    }
  }, error = function(e) NULL)
  
  # Method 2: Try looking in colData
  cd <- as.data.frame(colData(sce_obj))
  possible_col_pairs <- list(c("x", "y"), c("X", "Y"), c("center_x", "center_y"))
  for (cols in possible_col_pairs) {
      if (all(cols %in% colnames(cd))) {
          cat(paste("✓ Found spatial coordinates in", obj_name, "colData:", paste(cols, collapse=","), "\n"))
          coords <- cd[, cols]
          colnames(coords) <- c("x", "y")
          return(coords)
      }
  }
  
  # Method 3: Try looking in reducedDims
  if (!is.null(reducedDimNames(sce_obj))) {
      for (red_dim in reducedDimNames(sce_obj)) {
          if (grepl("spatial", tolower(red_dim))) {
              cat(paste("✓ Found spatial coordinates in", obj_name, "reducedDim:", red_dim, "\n"))
              coords <- as.data.frame(reducedDim(sce_obj, red_dim))
              colnames(coords)[1:2] <- c("x", "y")
              return(coords)
          }
      }
  }
  
  stop(paste("Error: No spatial coordinates found in", obj_name))
}

# Extract spatial coordinates
spatial_coords <- get_spatial_coords(dss9_data, "DSS9 Data")

# _________________
# Extract Metadata
# _________________

cat("\nExtracting metadata for annotations...\n")

metadata <- colData(dss9_data)

# Create plotting dataframe
plot_df <- data.frame(
    x = spatial_coords$x,
    y = spatial_coords$y,
    azimuth_annot = metadata$azimuth_celltype_l1_merged,
    azimuth_score = metadata$azimuth_mapping_score,
    singler_annot = metadata$singler_annotations_merged,
    singler_score = NA,
    celltypist_annot = metadata$celltypist_predicted_labels_merged,
    celltypist_score = metadata$celltypist_conf_score,
    scarches_annot = metadata$scvi_prediction_merged,
    scarches_score = metadata$scvi_confidence
)

# Extract SingleR scores (from matrix)
if (is.matrix(metadata$singler_scores)) {
    singler_scores_matrix <- metadata$singler_scores
    singler_assigned <- metadata$singler_annotations
    plot_df$singler_score <- sapply(1:nrow(metadata), function(i) {
        cell_type <- singler_assigned[i]
        singler_scores_matrix[i, cell_type]
    })
} else {
    plot_df$singler_score <- metadata$singler_scores
}

cat(sprintf("✓ Metadata extracted for %d cells.\n", nrow(plot_df)))

# ____________
# Setup Theme 
# ____________

# Text sizes massively increased
merfish_theme <- theme_void() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 26, color = "#1F2937", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.0, "cm"),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        aspect.ratio = 1
    )

# _______________
# Color Palettes
# _______________

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
    "Other"               = "#E8E8E8"
)

point_size <- 0.3
point_alpha <- 0.6

# ___________________
# Spatial Visualizations
# ___________________

cat("\nCreating spatial visualizations for annotation methods...\n")

# Base theme to elongate the continuous color bars on the right-side plots
continuous_leg_theme <- theme(
    legend.key.height = unit(2.0, "cm"),
    legend.key.width = unit(0.7, "cm")
)

# ________________
# Row 1: Azimuth
# ________________

cat("Creating Azimuth plots...\n")

p_s1a <- ggplot(plot_df, aes(x = x, y = y, color = azimuth_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    ggtitle("Azimuth - Cell Type") +
    merfish_theme

# Format the master cell type legend for the BOTTOM of the figure
p_legend_format <- p_s1a + 
    theme(
        legend.position = "bottom",
        legend.box.margin = margin(20, 0, 0, 0) # Add some breathing room above it
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 6), 
        nrow = 3,        # Spread horizontally across 3 rows
        byrow = TRUE
    ))

# Extract the master legend
shared_celltype_legend <- get_legend(p_legend_format)

p_s1b <- ggplot(plot_df, aes(x = x, y = y, color = azimuth_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Prediction Score", option = "magma") + 
    ggtitle("Azimuth - Prediction Score") +
    merfish_theme + continuous_leg_theme


# ________________
# Row 2: SingleR
# ________________

cat("Creating SingleR plots...\n")

p_s2a <- ggplot(plot_df, aes(x = x, y = y, color = singler_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    ggtitle("SingleR - Cell Type") +
    merfish_theme

p_s2b <- ggplot(plot_df, aes(x = x, y = y, color = singler_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Delta Score", option = "magma") + 
    ggtitle("SingleR - Delta Score") +
    merfish_theme + continuous_leg_theme

# ________________
# Row 3: CellTypist
# ________________

cat("Creating CellTypist plots...\n")

p_s3a <- ggplot(plot_df, aes(x = x, y = y, color = celltypist_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    ggtitle("CellTypist - Cell Type") +
    merfish_theme

p_s3b <- ggplot(plot_df, aes(x = x, y = y, color = celltypist_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Confidence Score", option = "magma") + 
    ggtitle("CellTypist - Confidence Score") +
    merfish_theme + continuous_leg_theme

# ________________
# Row 4: scArches
# ________________

cat("Creating scArches plots...\n")

p_s4a <- ggplot(plot_df, aes(x = x, y = y, color = scarches_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    ggtitle("scArches - Cell Type") +
    merfish_theme

p_s4b <- ggplot(plot_df, aes(x = x, y = y, color = scarches_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Uncertainty Score", option = "magma") + 
    ggtitle("scArches - Uncertainty Score") +
    merfish_theme + continuous_leg_theme

cat("✓ All spatial plots created.\n")

# _______________________
# Combine into Grid
# _______________________

cat("\nCombining plots...\n")

# Create the main 4x2 grid of plots, aggressively stripping all cell-type legends
main_plot_grid <- plot_grid(
    p_s1a + theme(legend.position = "none"), p_s1b,
    p_s2a + theme(legend.position = "none"), p_s2b,
    p_s3a + theme(legend.position = "none"), p_s3b,
    p_s4a + theme(legend.position = "none"), p_s4b,
    nrow = 4, ncol = 2, labels = "AUTO",
    label_size = 32, label_fontface = "bold",
    label_x = 0.02, label_y = 0.98,
    align = "v", axis = "lr",
    rel_widths = c(1, 1.25) # Slightly wider right column to fit the score legends naturally
)

# Stack the main grid ON TOP, and the shared legend ON THE BOTTOM
fig_merfish_s2 <- plot_grid(
    main_plot_grid, 
    shared_celltype_legend,
    nrow = 2, 
    rel_heights = c(1, 0.12) # 10-12% of the vertical space dedicated to the bottom legend
)

cat("✓ Combined figure created.\n")

# _________________
# Display and Save
# _________________

cat("\nDisplaying and saving plots...\n")
print(fig_merfish_s2)

dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)
# Scaled up total image size to safely accommodate the bigger text and the bottom legend
ggsave("figures/supp/merfish/Fig_S2_annotations_scores.png", fig_merfish_s2, width = 18, height = 26, dpi = 600)
cat("✓ Plot saved to figures/supp/merfish/\n")