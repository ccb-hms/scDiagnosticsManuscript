# ---------------------------------------------
# MERFISH Mouse Colon IBD - Spatial Figures
# ---------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scDiagnostics)
library(viridis)

# ---------------------------------------------

# _____________
# Data Loading
# _____________

cat("\n--- Loading Processed SCE Objects ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Setting the seed
set.seed(0)

# __________________________
# Define Metadata Mappings
# __________________________

cat("\n--- Mapping Cell Types ---\n")

if (!"tier2_merged" %in% colnames(colData(dss9_data))) {
  stop("Column 'tier2_merged' not found in dss9_data. Please run the annotation integration pipeline first.")
}

# Use the pre-calculated merged column for plotting categories
dss9_data$broad_panel_a <- dss9_data$tier2_merged
healthy_data$broad_panel_a <- healthy_data$tier2_merged

cat("✓ Using 'tier2_merged' as the source for cell type categories.\n")

# __________________
# Anomaly Detection 
# __________________

cat("\n--- Running Anomaly Detection ---\n")

# Define analysis class based on the merged columns
healthy_data$analysis_class <- ifelse(healthy_data$broad_panel_a == "Fibroblast", "Fibroblast_Lineage", "Ignore")
dss9_data$analysis_class <- ifelse(dss9_data$broad_panel_a %in% c("Fibroblast", "Inflamed Fibroblast"), "Fibroblast_Lineage", "Ignore")

anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "analysis_class", 
    ref_cell_type_col = "analysis_class", 
    cell_types = "Fibroblast_Lineage",
    pc_subset = 1:3, 
    anomaly_threshold = 0.5,
    max_cells_query = NULL,
    max_cells_ref = NULL
)

# Extract anomalous barcodes
anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# __________________________________
# Prepare Data & Calculate Metrics
# __________________________________

plot_df <- data.frame(
    barcode = colnames(dss9_data),
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    tier2 = dss9_data$tier2,
    broad_type = dss9_data$broad_panel_a # Using merged type
)

plot_df$anomaly_status <- ifelse(plot_df$barcode %in% anomalous_barcodes, "Anomalous", "Typical")

# Add ECM Score (Loss of Homeostasis Signature)
# Define standard fibrosis/ECM genes often found in MERFISH panels
target_ecm_genes <- c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Timp1", "Sparc", "Mmp2", "Vim")
available_genes <- intersect(target_ecm_genes, rownames(dss9_data))

cat(sprintf("Calculating ECM score using %d available genes: %s\n", 
            length(available_genes), paste(available_genes, collapse=", ")))

if(length(available_genes) > 0) {
    signature_expression <- assay(dss9_data, "logcounts")[available_genes, , drop = FALSE]
    
    # 1. Calculate raw mean
    raw_score <- colMeans(signature_expression)
    
    # 2. Cap at 1.5 (Winsorize) to prevent outliers washing out the plot
    capped_score <- pmin(raw_score, 1.5)
    
    # 3. Scale by 2 (Range is roughly 0 to 3 for plotting)
    plot_df$ecm_score <- capped_score * 2
    
} else {
    warning("No ECM genes found in dataset rows. ECM score will be 0.")
    plot_df$ecm_score <- 0
}

# METRICS REPORTING
cat("\n------------------------------------------------\n")
cat("      ANOMALY DETECTION PERFORMANCE METRICS       \n")
cat("      (Evaluated on Fibroblast Lineage Only)      \n")
cat("------------------------------------------------\n")

fibro_subset <- plot_df |> 
    filter(broad_type %in% c("Fibroblast", "Inflamed Fibroblast"))

TP <- sum(fibro_subset$broad_type == "Inflamed Fibroblast" & fibro_subset$anomaly_status == "Anomalous")
FN <- sum(fibro_subset$broad_type == "Inflamed Fibroblast" & fibro_subset$anomaly_status == "Typical")
TN <- sum(fibro_subset$broad_type == "Fibroblast" & fibro_subset$anomaly_status == "Typical")
FP <- sum(fibro_subset$broad_type == "Fibroblast" & fibro_subset$anomaly_status == "Anomalous")

sensitivity <- TP / (TP + FN) * 100
specificity <- TN / (TN + FP) * 100
accuracy <- (TP + TN) / (TP + TN + FP + FN) * 100

cat(sprintf("Total Inflamed Fibroblasts (Ground Truth): %d\n", (TP + FN)))
cat(sprintf("  - Identified as Anomalous (TP):          %d\n", TP))
cat(sprintf("  - Identified as Typical (FN):            %d\n", FN))
cat(sprintf("Total Normal Fibroblasts (Ground Truth):   %d\n", (TN + FP)))
cat(sprintf("  - Identified as Typical (TN):            %d\n", TN))
cat(sprintf("  - Identified as Anomalous (FP):          %d\n", FP))
cat("\n")
cat(sprintf("SENSITIVITY (Recall):  %.2f%%\n", sensitivity))
cat(sprintf("SPECIFICITY:           %.2f%%\n", specificity))
cat(sprintf("ACCURACY:              %.2f%%\n", accuracy))
cat("------------------------------------------------\n")

# _____________
# Visual Setup
# _____________

point_size <- 0.15 
bg_color <- "grey90"
common_theme <- theme_void() + 
    theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size=9, face="bold"),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.4, "cm")
    )

# --- PANEL A ---

# 1. Define Color Palette 
panel_a_colors <- c(
    # --- Fibroblasts ---
    "Inflamed Fibroblast" = "#A50F15",  
    "Fibroblast"          = "#FB9A99",  
    
    # --- Smooth Muscle ---
    "Inflamed SMC"        = "#54278F",  
    "Smooth Muscle"       = "#BC80BD",  
    
    # --- Epithelial ---
    "Inflamed Epithelial" = "#E7298A",  
    "Epithelial"          = "#F1B6DA",  
    "Stem/TA"             = "#008B8B",  
    
    # --- Immune Cells ---
    "Other Immune"        = "#FDB462",  
    "Neutrophil"          = "#4682B4",  
    
    # --- Others ---
    "Pericyte/FRC"        = "#8C510A",  
    "Endothelial"         = "#B2DF8A",  
    "Enteric Nervous"     = "#FFFF99",  
    "ICC"                 = "#8DD3C7",  
    "Adipose"             = "#E5C494",  
    "Other"               = bg_color
)

# 2. Define Legend Order 
legend_order <- c(
    "Fibroblast", "Inflamed Fibroblast",
    "Smooth Muscle", "Inflamed SMC",
    "Epithelial", "Stem/TA", "Inflamed Epithelial", # Grouped together
    "Endothelial", "Pericyte/FRC", 
    "Neutrophil", "Other Immune", 
    "Enteric Nervous", "ICC", "Adipose", "Other"
)

# 3. Define Z-Order (Stem/TA goes in Healthy Layer)
z_order_levels <- c(
    "Other", "Adipose", "Enteric Nervous", "ICC", "Endothelial", "Pericyte/FRC",
    "Other Immune", "Epithelial", "Stem/TA", "Smooth Muscle", "Fibroblast", # Healthy Layer
    "Neutrophil", "Inflamed Epithelial", "Inflamed SMC", "Inflamed Fibroblast" # Top Layer
)

# Apply Factor for Legend
plot_df$broad_type <- factor(plot_df$broad_type, levels = legend_order)

# Create Sorted Dataframe for Plotting (Z-Order)
plot_df_a <- plot_df |>
    mutate(z_sort = factor(broad_type, levels = z_order_levels)) |>
    arrange(z_sort)

panel_a <- ggplot(plot_df_a, aes(x=x, y=y, color=broad_type)) +
    geom_point(size = point_size) +
    scale_color_manual(values = panel_a_colors, name = "Cell Type") +
    labs(title = "A. Spatial Cell Type Map") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3), ncol=4))

# --- PANEL B: Fibroblast State (Annotation) ---

# Matching Panel A Colors
panel_b_df <- plot_df
panel_b_df$display_group <- "Other"
panel_b_df$display_group[panel_b_df$broad_type == "Fibroblast"] <- "Healthy Fibroblast"
panel_b_df$display_group[panel_b_df$broad_type == "Inflamed Fibroblast"] <- "Inflamed Fibroblast (IAF)"

panel_b_df$display_group <- factor(panel_b_df$display_group, 
                                   levels = c("Other", "Healthy Fibroblast", "Inflamed Fibroblast (IAF)"))
panel_b_df <- panel_b_df |> arrange(display_group)

panel_b_colors <- c(
    "Healthy Fibroblast"        = "#FB9A99", 
    "Inflamed Fibroblast (IAF)" = "#A50F15", 
    "Other"                     = bg_color
)

panel_b <- ggplot(panel_b_df, aes(x=x, y=y, color=display_group)) +
    geom_point(size = point_size) +
    scale_color_manual(values = panel_b_colors, name = "Fibroblast State") +
    labs(title = "B. Fibroblast State (Annotation)") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3)))

# --- PANEL C: Anomalies ---

panel_c_df <- plot_df
panel_c_df$display_group <- "Other"
is_fibro_lineage <- panel_c_df$broad_type %in% c("Fibroblast", "Inflamed Fibroblast")
panel_c_df$display_group[is_fibro_lineage & panel_c_df$anomaly_status == "Typical"] <- "Typical"
panel_c_df$display_group[is_fibro_lineage & panel_c_df$anomaly_status == "Anomalous"] <- "Anomalous"

panel_c_df$display_group <- factor(panel_c_df$display_group, 
                                   levels = c("Other", "Typical", "Anomalous"))
panel_c_df <- panel_c_df |> arrange(display_group)

panel_c_colors <- c(
    "Typical"   = "#33A02C", 
    "Anomalous" = "#E31A1C", 
    "Other"     = bg_color
)

panel_c <- ggplot(panel_c_df, aes(x=x, y=y, color=display_group)) +
    geom_point(size = point_size) +
    scale_color_manual(values = panel_c_colors, name = "scDiagnostics Prediction") +
    labs(title = "C. Detected Anomalies") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3)))

# --- PANEL D: ECM Score ---

GLOBAL_ECM_LIMIT <- 3.0

panel_d_df <- plot_df |> arrange(ecm_score)

panel_d <- ggplot() +
    geom_point(data = subset(panel_d_df, !broad_type %in% c("Fibroblast", "Inflamed Fibroblast")), 
               aes(x=x, y=y), color = bg_color, size = point_size) +
    geom_point(data = subset(panel_d_df, broad_type %in% c("Fibroblast", "Inflamed Fibroblast")), 
               aes(x=x, y=y, color = ecm_score), size = point_size) +
    scale_color_viridis_c(
        option = "plasma", 
        direction = 1, 
        name = "ECM Score",
        limits = c(0, GLOBAL_ECM_LIMIT),
        oob = scales::squish
    ) +
    labs(title = "D. Loss of Homeostasis Signature") +
    common_theme

# ________________
# Assemble & Save
# ________________

cat("\n--- Assembling Final Figure ---\n")
final_figure <- (panel_a + panel_b) / (panel_c + panel_d) + 
    plot_annotation(
        title = "Spatial Validation of Anomalous Fibroblasts (MERFISH)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

print(final_figure)

dir.create("figures/merfish", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/merfish/spatial_plots.png", 
       plot = final_figure, width = 14, height = 14, units = "in", dpi = 600)

cat("\n✓ Spatial plot saved.\n")