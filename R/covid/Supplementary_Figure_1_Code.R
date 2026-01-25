# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 1 Code
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

theme_set(theme_minimal() + theme(
  axis.title = element_text(size = 11, face = "bold"),
  axis.text = element_text(size = 10),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 9, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 10, face = "bold")
))

# Extract UMAP coordinates for all cells
umap_coords <- reducedDim(covid_data, "UMAP_scVI")

# _________________________________________
# Figure S1A: Azimuth - Cell Type (Merged)
# _________________________________________

umap_data_azimuth_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$azimuth_celltype_l1_merged
)

fig_s1a <- ggplot(umap_data_azimuth_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    guides(color = guide_legend(ncol = 1)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Azimuth - Cell Type") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _______________________________________
# Figure S1B: Azimuth - Prediction Score
# _______________________________________

umap_data_azimuth_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    PredictionScore = covid_data$azimuth_mapping_score
)

fig_s1b <- ggplot(umap_data_azimuth_score, aes(x = UMAP1, y = UMAP2, color = PredictionScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Prediction Score") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Azimuth - Prediction Score") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _________________________________________
# Figure S1C: SingleR - Cell Type (Merged)
# _________________________________________

umap_data_singler_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$singler_annotations_merged
)

fig_s1c <- ggplot(umap_data_singler_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    guides(color = guide_legend(ncol = 1)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("SingleR - Cell Type") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# ________________________________________
# Figure S1D: SingleR - Delta Score
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

fig_s1d <- ggplot(umap_data_singler_score, aes(x = UMAP1, y = UMAP2, color = DeltaScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Delta Score") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("SingleR - Delta Score") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# ____________________________________________
# Figure S1E: CellTypist - Cell Type (Merged)
# ____________________________________________

umap_data_celltypist_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$celltypist_predicted_labels_merged
)

fig_s1e <- ggplot(umap_data_celltypist_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    guides(color = guide_legend(ncol = 1)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("CellTypist - Cell Type") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _________________________________________
# Figure S1F: CellTypist - Confidence Score
# _________________________________________

umap_data_celltypist_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    ConfidenceScore = covid_data$celltypist_conf_score
)

fig_s1f <- ggplot(umap_data_celltypist_score, aes(x = UMAP1, y = UMAP2, color = ConfidenceScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Confidence Score") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("CellTypist - Confidence Score") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _________________________________________
# Figure S1G: scArches - Cell Type (Merged)
# _________________________________________

umap_data_scarches_anno <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = covid_data$scvi_prediction_merged
)

fig_s1g <- ggplot(umap_data_scarches_anno, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(values = cell_type_colors, na.value = "grey50") +
    guides(color = guide_legend(ncol = 1)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("scArches - Cell Type") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _________________________________________
# Figure S1H: scArches - Uncertainty Score
# _________________________________________

umap_data_scarches_score <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    UncertaintyScore = covid_data$scvi_confidence
)

fig_s1h <- ggplot(umap_data_scarches_score, aes(x = UMAP1, y = UMAP2, color = UncertaintyScore)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient(low = "#5B7C99", high = "#C1666B", name = "Uncertainty Score") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("scArches - Uncertainty Score") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 9, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 8, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(2, 0.5, 2, 0.5, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

# _________________
# Combine 4x2 grid
# _________________

fig_s1_combined <- plot_grid(fig_s1a, fig_s1b, 
                             fig_s1c, fig_s1d,
                             fig_s1e, fig_s1f,
                             fig_s1g, fig_s1h,
                             nrow = 4, ncol = 2, labels = "AUTO", 
                             label_size = 11, label_fontface = "bold",
                             label_x = 0.02, label_y = 0.98,
                             rel_widths = c(0.9, 0.9))

ggsave("figures/supp/covid/Fig_S1_umaps.png", fig_s1_combined, width = 15, height = 16, dpi = 600)
print("Figure S1 saved")

# ________
# Summary
# ________

print("Supplementary Figure S1 (All Cells) complete!")
print("Saved in: figures/supp/")