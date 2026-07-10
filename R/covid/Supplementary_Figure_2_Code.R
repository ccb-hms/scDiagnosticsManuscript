# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 2 Code
# -----------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(dplyr)

# -----------------------------------------------

# __________
# Load data
# __________

covid_data <- readRDS("data/covid/covid_data_sce.rds")

# _________________________
# Cell Type Color Palette
# _________________________

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

# ____________
# Setup Theme
# ____________

theme_set(theme_minimal(base_size = 20) + theme(
  axis.title = element_text(size = 20, face = "bold"),
  axis.text = element_text(size = 16),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 20, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 18, face = "bold")
))

# Extract UMAP coordinates for all cells
umap_coords <- reducedDim(covid_data, "UMAP_scVI")

# _____________________________________________
# Custom Themes for Left & Right Columns
# _____________________________________________

# LEFT: No legend, big top margin, reduced right margin to pull right column closer
theme_left <- theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "#1F2937", family = "sans"),
    axis.title = element_text(size = 20, face = "bold", color = "#374151", family = "sans"),
    axis.text = element_text(size = 16, color = "#4B5563", family = "sans"),
    legend.position = "none", 
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), # Removed the frame box
    axis.line = element_line(color = "#374151", linewidth = 0.5), # Added back just the X/Y axes lines
    plot.margin = margin(t = 35, r = 5, b = 10, l = 10, "pt"), 
    
    # --- CRITICAL FIX: Force pure white backgrounds ---
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    aspect.ratio = 1
)

# RIGHT: Keeps score legend, big top margin, reduced left margin to sit flush with left column
theme_right <- theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "#1F2937", family = "sans"),
    axis.title = element_text(size = 20, face = "bold", color = "#374151", family = "sans"),
    axis.text = element_text(size = 16, color = "#4B5563", family = "sans"),
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.6, "cm"),
    legend.background = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), # Removed the frame box
    axis.line = element_line(color = "#374151", linewidth = 0.5), # Added back just the X/Y axes lines
    plot.margin = margin(t = 35, r = 10, b = 10, l = 0, "pt"), 
    
    # --- CRITICAL FIX: Force pure white backgrounds ---
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    aspect.ratio = 1
)


# _________________________________________
# Figure S2A: Azimuth - Cell Type (Merged)
# _________________________________________

umap_data_azimuth_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$azimuth_celltype_l1_merged
)

fig_s2a <- ggplot(umap_data_azimuth_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Azimuth - Cell Type") + theme_left

# _______________________________________
# Figure S2B: Azimuth - Prediction Score
# _______________________________________

umap_data_azimuth_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    PredictionScore = covid_data$azimuth_mapping_score
)

fig_s2b <- ggplot(umap_data_azimuth_score, aes(x = UMAP1, y = UMAP2, color = PredictionScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Prediction Score") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("Azimuth - Prediction Score") + theme_right

# _________________________________________
# Figure S2C: SingleR - Cell Type (Merged)
# _________________________________________

umap_data_singler_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$singler_annotations_merged
)

fig_s2c <- ggplot(umap_data_singler_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("SingleR - Cell Type") + theme_left

# ________________________________________
# Figure S2D: SingleR - Delta Score
# ________________________________________

singler_scores_matrix <- covid_data$singler_scores
singler_assigned <- covid_data$singler_annotations

singler_conf <- sapply(1:ncol(covid_data), function(i) {
    cell_type <- singler_assigned[i]
    singler_scores_matrix[i, cell_type]
})

umap_data_singler_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    DeltaScore = singler_conf
)

fig_s2d <- ggplot(umap_data_singler_score, aes(x = UMAP1, y = UMAP2, color = DeltaScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Delta Score") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("SingleR - Delta Score") + theme_right

# ____________________________________________
# Figure S2E: CellTypist - Cell Type (Merged)
# ____________________________________________

umap_data_celltypist_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$celltypist_predicted_labels_merged
)

fig_s2e <- ggplot(umap_data_celltypist_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("CellTypist - Cell Type") + theme_left

# _________________________________________
# Figure S2F: CellTypist - Confidence Score
# _________________________________________

umap_data_celltypist_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    ConfidenceScore = covid_data$celltypist_conf_score
)

fig_s2f <- ggplot(umap_data_celltypist_score, aes(x = UMAP1, y = UMAP2, color = ConfidenceScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Confidence Score") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("CellTypist - Confidence Score") + theme_right

# _________________________________________
# Figure S2G: scArches - Cell Type (Merged)
# _________________________________________

umap_data_scarches_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$scvi_prediction_merged
)

fig_s2g <- ggplot(umap_data_scarches_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("scArches - Cell Type") + theme_left

# _________________________________________
# Figure S2H: scArches - Uncertainty Score
# _________________________________________

umap_data_scarches_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    UncertaintyScore = covid_data$scvi_confidence
)

fig_s2h <- ggplot(umap_data_scarches_score, aes(x = UMAP1, y = UMAP2, color = UncertaintyScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Uncertainty Score") +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("scArches - Uncertainty Score") + theme_right

# _________________
# Combine 4x2 grid
# _________________

# 1. Extract the shared Cell Type legend
dummy_plot <- ggplot(umap_data_azimuth_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 5) + scale_color_manual(values = cell_type_colors) +
    theme_void() + 
    theme(
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 18),
        # Ensure the legend dummy plot doesn't have transparency issues either
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(nrow = 4, override.aes = list(size = 6)))

shared_legend <- cowplot::get_legend(dummy_plot)

# 2. Build the main grid 
grid_plots <- plot_grid(
    fig_s2a, fig_s2b, 
    fig_s2c, fig_s2d,
    fig_s2e, fig_s2f,
    fig_s2g, fig_s2h,
    nrow = 4, ncol = 2, labels = "AUTO", 
    label_size = 26, label_fontface = "bold",
    label_x = 0, label_y = 1, # Snaps labels into the top margin we created
    rel_widths = c(1, 1.15)
)

# 3. Append the shared legend at the bottom
fig_s2_combined <- plot_grid(grid_plots, shared_legend, ncol = 1, rel_heights = c(1, 0.12))

# *** CRITICAL FIX: bg = "white" prevents transparent PNG artifacts ***
ggsave("figures/supp/covid/Fig_S2_umaps.png", fig_s2_combined, width = 15, height = 24, dpi = 600, bg = "white")
print("Figure S2 saved")

# ________
# Summary
# ________

print("Supplementary Figure S2 (All Cells) complete!")
print("Saved in: figures/supp/")