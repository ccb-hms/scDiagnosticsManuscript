# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 1 Code 
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

# _____________________________________
# Create Combined Plotting Data Frame
# _____________________________________

cat("\nCombining data for plotting...\n")

if (!"tier2_merged" %in% colnames(colData(normal_data)) || !"tier2_merged" %in% colnames(colData(dss9_data))) {
  stop("Column 'tier2_merged' missing. Run the merging script first.")
}

meta_normal <- as.data.frame(colData(normal_data)) |>
    select(tier2_merged) |>
    mutate(disease_status = "Healthy",
           cell_type = tier2_merged,
           barcode = rownames(colData(normal_data)))

meta_dss9 <- as.data.frame(colData(dss9_data)) |>
    select(tier2_merged) |>
    mutate(disease_status = "DSS9 (Colitis)",
           cell_type = tier2_merged,
           barcode = rownames(colData(dss9_data)))

plot_df <- bind_rows(
    bind_cols(umap_normal, meta_normal),
    bind_cols(umap_dss9, meta_dss9)
)
cat(sprintf("✓ Combined plotting data frame created with %d cells.\n", nrow(plot_df)))

# _______________________
# ECM Score Computation
# _______________________

cat("\nComputing ECM Homeostasis Score...\n")

ecm_genes <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")

common_genes <- intersect(rownames(normal_data), rownames(dss9_data))
available_ecm_genes <- intersect(ecm_genes, common_genes)

cat("Requested genes:", paste(ecm_genes, collapse=", "), "\n")
cat("Available genes:", paste(available_ecm_genes, collapse=", "), "\n")

if(length(available_ecm_genes) == 0) {
    warning("None of the ECM genes were found in the dataset! Setting score to 0.")
    plot_df$ECM_Score <- 0
} else {
    expr_normal <- assay(normal_data, "logcounts")[available_ecm_genes, , drop=FALSE]
    expr_dss9   <- assay(dss9_data, "logcounts")[available_ecm_genes, , drop=FALSE]
    
    scores_normal <- colMeans(as.matrix(expr_normal), na.rm = TRUE)
    scores_dss9   <- colMeans(as.matrix(expr_dss9), na.rm = TRUE)
    
    all_scores <- c(scores_normal, scores_dss9)
    capped_scores <- pmin(all_scores, 1.5)
    final_scores <- capped_scores * 2
    
    plot_df$ECM_Score <- final_scores[match(plot_df$barcode, names(final_scores))]
    
    cat("✓ ECM scores computed (Range: ", min(plot_df$ECM_Score), " - ", max(plot_df$ECM_Score), ")\n")
}

# ____________
# Setup Theme 
# ____________

merfish_theme <- theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 26, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 22, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 18, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.0, "cm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),                                # No box
        axis.line = element_line(color = "#374151", linewidth = 0.5),  # Just L-shaped axes
        plot.margin = margin(20, 20, 20, 20, "pt"),
        panel.background = element_rect(fill = "white", color = NA),   # Force pure white inner
        plot.background = element_rect(fill = "white", color = NA),    # Force pure white outer margins!
        aspect.ratio = 1
    )

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

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

disease_colors <- c("Healthy" = "#5DADE2", "DSS9 (Colitis)" = "#FF6347")

point_size <- 0.5
point_alpha <- 0.6

# Plot 1: Disease Status
cat("Creating UMAP plot by disease status...\n")  
p_disease <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = disease_status)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = disease_colors, name = "Condition") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Disease Status") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 5))) 

# Plot 2: Cell Type
cat("Creating UMAP plot by cell type...\n")
types_to_keep <- plot_df |> count(cell_type) |> filter(n >= 10) |> pull(cell_type)
plot_df_filtered <- plot_df |> filter(cell_type %in% types_to_keep)

p_celltype <- ggplot(plot_df_filtered, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = merfish_colors, name = "Cell Type", na.value = "grey80") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Cell Type") +
    merfish_theme +
    guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) 

# Plot 3: ECM Score
cat("Creating UMAP plot by ECM Score...\n")
p_ecm <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = ECM_Score)) +
    geom_point(size = point_size, alpha = 0.8) + 
    scale_color_viridis_c(name = "ECM Score", option = "magma") + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("ECM Homeostasis Score") +
    merfish_theme +
    theme(
        legend.key.height = unit(1.5, "cm"), 
        legend.key.width = unit(0.5, "cm")
    )

cat("✓ UMAP plots created.\n")

# _______________________
# Combine into 3x1 Grid
# _______________________

cat("\nCombining plots into 3x1 figure...\n")

fig_merfish_combined <- plot_grid(
    p_disease, p_celltype, p_ecm,
    nrow = 3, ncol = 1, labels = "AUTO", 
    label_size = 32, label_fontface = "bold", 
    label_x = 0.02, label_y = 0.98,
    align = "v", axis = "lr" 
)

cat("✓ Combined figure created.\n")

# _________________
# Display and Save
# _________________

cat("\nDisplaying and saving plots...\n")
print(fig_merfish_combined)

dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)

# The crucial fix here: bg = "white" prevents transparent PNG artifacts
ggsave("figures/supp/merfish/Fig_S1_umaps_exploratory.png", fig_merfish_combined, width = 12, height = 30, dpi = 600, bg = "white")
cat("✓ Plots saved to figures/supp/merfish/\n")