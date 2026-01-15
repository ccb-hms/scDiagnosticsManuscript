# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 2 Code
# -------------------------------------------------------

library(SingleCellExperiment)
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

# Extract UMAP coordinates
umap_coords <- get_umap_coords(dss9_data, "DSS9 Data")

# _________________
# Extract Metadata
# _________________

cat("\nExtracting metadata for annotations...\n")

metadata <- colData(dss9_data)

# Create plotting dataframe
plot_df <- data.frame(
    UMAP1 = umap_coords$UMAP1,
    UMAP2 = umap_coords$UMAP2,
    singler_annot = metadata$singler_annotations_merged,
    singler_score = NA,
    azimuth_annot = metadata$azimuth_celltype_l1_merged,
    azimuth_score = metadata$azimuth_mapping_score,
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

merfish_theme <- theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 2, 2, 2, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
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
# UMAP Visualizations
# ___________________

cat("\nCreating UMAP visualizations for annotation methods...\n")

# ________________
# Row 1: SingleR
# ________________

cat("Creating SingleR plots...\n")

p_s1a <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = singler_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("SingleR - Cell Type") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))

p_s1b <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = singler_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Delta Score", option = "magma") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("SingleR - Delta Score") +
    merfish_theme

# ________________
# Row 2: Azimuth
# ________________

cat("Creating Azimuth plots...\n")

p_s2a <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = azimuth_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Azimuth - Cell Type") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))

p_s2b <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = azimuth_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Prediction Score", option = "magma") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Azimuth - Prediction Score") +
    merfish_theme

# ________________
# Row 3: CellTypist
# ________________

cat("Creating CellTypist plots...\n")

p_s3a <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = celltypist_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("CellTypist - Cell Type") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))

p_s3b <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = celltypist_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Confidence Score", option = "magma") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("CellTypist - Confidence Score") +
    merfish_theme

# ________________
# Row 4: scArches
# ________________

cat("Creating scArches plots...\n")

p_s4a <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = scarches_annot)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("scArches - Cell Type") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))

p_s4b <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = scarches_score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "Uncertainty Score", option = "magma") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("scArches - Uncertainty Score") +
    merfish_theme

cat("✓ All UMAP plots created.\n")

# _______________________
# Combine into 4x2 Grid
# _______________________

cat("\nCombining plots into 4x2 figure...\n")

fig_merfish_s2 <- plot_grid(
    p_s1a, p_s1b,
    p_s2a, p_s2b,
    p_s3a, p_s3b,
    p_s4a, p_s4b,
    nrow = 4, ncol = 2, labels = "AUTO",
    label_size = 12, label_fontface = "bold",
    label_x = 0.02, label_y = 0.98,
    rel_widths = c(1.1, 1)
)

cat("✓ Combined figure created.\n")

# _________________
# Display and Save
# _________________

cat("\nDisplaying and saving plots...\n")
print(fig_merfish_s2)

dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/supp/merfish/Fig_S2_annotations_scores.png", fig_merfish_s2, width = 14, height = 16, dpi = 600)
cat("✓ Plot saved to figures/supp/merfish/\n")
