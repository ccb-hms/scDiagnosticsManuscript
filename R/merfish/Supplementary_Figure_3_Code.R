# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 3 Code
# -------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
library(scDiagnostics)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\nLoading datasets for inflamed fibroblast analysis...\n")
normal_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data   <- readRDS("data/merfish/dss9_data.rds")

cat("Datasets loaded.\n")

# ______
# Setup
# ______

inflamed_fb_mask <- dss9_data$tier2_merged == "Inflamed Fibroblast"
n_inflamed_fb <- sum(inflamed_fb_mask)
cat(sprintf("\n✓ Total inflamed fibroblasts: %d\n", n_inflamed_fb))

all_cell_types <- unique(normal_data$tier2_merged)
cat(sprintf("✓ Cell types in reference: %s\n", paste(all_cell_types, collapse = ", ")))

inflamed_fb_cellnames <- colnames(dss9_data)[inflamed_fb_mask]

# ______________________________________
# Run Anomaly Detection for All Methods
# ______________________________________

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("ANOMALY DETECTION\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n")

cat("\nDetecting anomalies for Azimuth...\n")
anomaly_azimuth <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "azimuth_celltype_l1_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for SingleR...\n")
anomaly_singler <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "singler_annotations_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for CellTypist...\n")
anomaly_celltypist <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "celltypist_predicted_labels_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for scArches...\n")
anomaly_scarches <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "scvi_prediction_merged",
    cell_types = all_cell_types[1:9],
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\n✓ All anomaly detections complete\n")

# _________________________________
# Create Plotting DataFrame
# _________________________________

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("PROCESSING RESULTS FOR HEATMAP\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n\n")

# Identify all unique predicted types seen across the 4 methods
all_pred_types <- sort(unique(c(
    dss9_data$azimuth_celltype_l1_merged[inflamed_fb_mask],
    dss9_data$singler_annotations_merged[inflamed_fb_mask],
    dss9_data$celltypist_predicted_labels_merged[inflamed_fb_mask],
    dss9_data$scvi_prediction_merged[inflamed_fb_mask]
)))

df_list <- list()

methods <- c("Azimuth", "SingleR", "CellTypist", "scArches")
pred_cols <- c("azimuth_celltype_l1_merged", "singler_annotations_merged", 
               "celltypist_predicted_labels_merged", "scvi_prediction_merged")
anomaly_results <- list(anomaly_azimuth, anomaly_singler, anomaly_celltypist, anomaly_scarches)

for (i in seq_along(methods)) {
    method_name <- methods[i]
    pred_col <- pred_cols[i]
    anomaly_res <- anomaly_results[[i]]
    
    names(dss9_data[[pred_col]]) <- colnames(dss9_data)
    pred_inflamed <- dss9_data[[pred_col]][inflamed_fb_mask]
    
    for (pred_type in all_pred_types) {
        mask_inflamed <- pred_inflamed == pred_type
        n_total <- sum(mask_inflamed, na.rm = TRUE)
        n_anom <- 0
        
        if (n_total > 0 && pred_type %in% names(anomaly_res)) {
            inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
            anomaly_for_type <- anomaly_res[[pred_type]][["query_anomaly"]]
            names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
            inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
            inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
            n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        }
        
        df_list[[length(df_list) + 1]] <- data.frame(
            Method = method_name,
            CellType = pred_type,
            Total = n_total,
            Anomalous = n_anom,
            stringsAsFactors = FALSE
        )
    }
}

plot_df <- do.call(rbind, df_list)

# Calculate Percentage Anomalous
plot_df$Pct <- ifelse(plot_df$Total > 0, (plot_df$Anomalous / plot_df$Total) * 100, NA)

# Create fraction text label (stacked vertically)
plot_df$Label <- ifelse(plot_df$Total > 0, sprintf("%d\n—\n%d", plot_df$Anomalous, plot_df$Total), "")

# Dynamic text color (Black text on bright yellow backgrounds, White text on dark purple backgrounds)
plot_df$TextColor <- ifelse(is.na(plot_df$Pct), "transparent", 
                            ifelse(plot_df$Pct > 65, "black", "white"))

# Lock factors so order is top-to-bottom on Y axis
plot_df$Method <- factor(plot_df$Method, levels = rev(c("Azimuth", "SingleR", "CellTypist", "scArches")))

cat("✓ Plotting dataframe ready.\n")

# ________________
# Create Plot
# ________________

cat("\nCreating combined heatmap...\n")

p_heatmap <- ggplot(plot_df, aes(x = CellType, y = Method)) +
    # The Heatmap tiles
    geom_tile(aes(fill = Pct), color = "white", linewidth = 1.5) +
    # The text overlay
    geom_text(aes(label = Label, color = TextColor), size = 6.5, fontface = "bold", lineheight = 0.8) +
    # Color scales
    scale_fill_viridis_c(
        option = "viridis", 
        na.value = "#F3F4F6", # Light grey for 0/0 cells so they visually disappear
        name = "Detection Rate\n(% Anomalous)", 
        limits = c(0, 100),
        labels = function(x) paste0(x, "%")
    ) +
    scale_color_identity() + # Uses the dynamic colors defined in the dataframe
    coord_fixed() + # Forces squares
    labs(
        title = "scDiagnostics Anomaly Detection on Inflamed Fibroblasts",
        subtitle = "Color = % of assigned cells flagged as anomalous  |  Text = Anomalous Flagged / Total Assigned",
        x = "Assigned Cell Type by Annotation Tool",
        y = ""
    ) +
    theme_minimal(base_size = 18) +
    theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0.5, color = "#1F2937", margin = margin(b = 10)),
        plot.subtitle = element_text(size = 20, hjust = 0.5, color = "#4B5563", margin = margin(b = 20)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black", face = "bold"),
        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
        axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 20)),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.height = unit(1.5, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 20, 20, "pt")
    )

dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/supp/merfish/Fig_S3_annotation_heatmap.png", p_heatmap, width = 20, height = 10, dpi = 600, bg = "white")

cat("\n✓ Combined heatmap saved to figures/supp/merfish/Fig_S3_annotation_heatmap.png\n")