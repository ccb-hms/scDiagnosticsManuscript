# -----------------------------------------------------------
# MERFISH Mouse Colon IBD - Final 4-Panel Spatial Figure
# Professional Colors & Detailed Metrics
# -----------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scDiagnostics)
library(viridis)

# -----------------------------------------------------------
# 1. Data Loading
# -----------------------------------------------------------
cat("\n--- Loading Processed SCE Objects ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Ensure annotations are characters
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)

set.seed(0)

# -----------------------------------------------------------
# 2. Define Metadata Mappings & Legend Order
# -----------------------------------------------------------
cat("\n--- Mapping Cell Types ---\n")

map_cell_types <- function(labels) {
  case_when(
    # --- FIBROBLAST (Blue Pair) ---
    grepl("^IAF", labels) ~ "Inflamed Fibroblast",
    grepl("^Fibro|^FRC|^Pericyte", labels) ~ "Fibroblast",
    
    # --- SMOOTH MUSCLE (Purple Pair) ---
    grepl("^IASMC", labels) ~ "Inflamed SMC",
    grepl("SMC", labels) ~ "Smooth Muscle",
    
    # --- EPITHELIAL (Rose/Berry Pair) ---
    # Corrected: IAE is grouped here
    grepl("^IAE", labels) ~ "Inflamed Epithelial",
    grepl("Colonocytes|Goblet|Stem|TA|Epithelial|M cells|Repair associated", labels) ~ "Epithelial",
    
    # --- OTHERS (Singletons) ---
    grepl("Endothelial|EC", labels) ~ "Endothelial",
    grepl("Glia|Neuron", labels) ~ "Enteric Nervous",
    grepl("ICC", labels) ~ "ICC",
    grepl("Adipose", labels) ~ "Adipose",
    grepl("B cell|T|Macrophage|Monocyte|Neutrophil|DC|ILC2|Plasma|Mast", labels) ~ "Immune",
    
    TRUE ~ "Other"
  )
}

dss9_data$broad_panel_a <- map_cell_types(dss9_data$tier2)
healthy_data$broad_panel_a <- map_cell_types(healthy_data$tier2)

# -----------------------------------------------------------
# 3. Anomaly Detection (Fibroblast Lineage Only)
# -----------------------------------------------------------
cat("\n--- Running Anomaly Detection ---\n")

# Isolate Fibroblasts for analysis
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

anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# -----------------------------------------------------------
# 4. Prepare Data & Calculate Metrics
# -----------------------------------------------------------
plot_df <- data.frame(
    barcode = colnames(dss9_data),
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    tier2 = dss9_data$tier2,
    broad_type = dss9_data$broad_panel_a
)

plot_df$anomaly_status <- ifelse(plot_df$barcode %in% anomalous_barcodes, "Anomalous", "Typical")

# Add ECM Score
healthy_fibro_signature_genes <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")
available_genes <- intersect(healthy_fibro_signature_genes, rownames(dss9_data))
if(length(available_genes) > 0) {
    signature_expression <- assay(dss9_data, "logcounts")[available_genes, , drop = FALSE]
    plot_df$ecm_score <- colMeans(signature_expression)
} else {
    plot_df$ecm_score <- 0
}

# --- METRICS REPORTING ---
cat("\n------------------------------------------------\n")
cat("      ANOMALY DETECTION PERFORMANCE METRICS       \n")
cat("      (Evaluated on Fibroblast Lineage Only)      \n")
cat("------------------------------------------------\n")

fibro_subset <- plot_df %>% 
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

# --- Enforce Legend Order ---
legend_order <- c(
    "Fibroblast", "Inflamed Fibroblast",
    "Smooth Muscle", "Inflamed SMC",
    "Epithelial", "Inflamed Epithelial",
    "Endothelial", "Immune", "Enteric Nervous", "ICC", "Adipose", "Other"
)
plot_df$broad_type <- factor(plot_df$broad_type, levels = legend_order)

# -----------------------------------------------------------
# 5. Generate Panels
# -----------------------------------------------------------
point_size <- 0.15 
bg_color <- "grey96"
common_theme <- theme_void() + 
    theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size=9, face="bold"),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.4, "cm")
    )

# --- PANEL A: PROFESSIONAL PALETTE ---
# --- PANEL A: UPDATED PROFESSIONAL PALETTE ---
panel_a_colors <- c(
    # --- PAIRED GROUPS (Fixed) ---
    # Fibroblasts (Blue)
    "Fibroblast" = "#A6CEE3",           # Light Blue
    "Inflamed Fibroblast" = "#1F78B4",  # Dark Blue
    
    # Smooth Muscle (Purple)
    "Smooth Muscle" = "#BC80BD",        # Light Lavender
    "Inflamed SMC" = "#7570B3",         # Deep Purple
    
    # Epithelial (Rose/Berry)
    "Epithelial" = "#FCCDE5",           # Pale Rose
    "Inflamed Epithelial" = "#C51B7D",  # Deep Magenta/Berry
    
    # --- SINGLETONS (Updated) ---
    # Immune: Light Orange/Gold (High contrast vs Blue, distinct from Pink)
    "Immune" = "#B2DF8A",               
    
    # Endothelial: Teal/Mint
    "Endothelial" = "#FDB462",          
    
    # Adipose: Tan (No Gray)
    "Adipose" = "#E5C494",              
    
    # Others
    "Enteric Nervous" = "#FFFF99",      # Pale Yellow
    "ICC" = "#8DD3C7",                  # Light Olive
    "Other" = bg_color
)

panel_a <- ggplot(plot_df, aes(x=x, y=y, color=broad_type)) +
    geom_point(size = point_size) +
    scale_color_manual(values = panel_a_colors, name = "Cell Type") +
    labs(title = "A. Spatial Cell Type Map") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3), ncol=4))

# --- PANEL B: Fibroblast State (Matches Panel A Blues) ---
panel_b_df <- plot_df
panel_b_df$display_group <- "Other"
panel_b_df$display_group[panel_b_df$broad_type == "Fibroblast"] <- "Healthy Fibroblast"
panel_b_df$display_group[panel_b_df$broad_type == "Inflamed Fibroblast"] <- "Inflamed Fibroblast (IAF)"
panel_b_df$display_group <- factor(panel_b_df$display_group, levels = c("Inflamed Fibroblast (IAF)", "Healthy Fibroblast", "Other"))

panel_b_colors <- c(
    "Healthy Fibroblast" = "#A6CEE3",        
    "Inflamed Fibroblast (IAF)" = "#1F78B4", 
    "Other" = bg_color
)

panel_b <- ggplot() +
    geom_point(data = subset(panel_b_df, display_group == "Other"), 
               aes(x=x, y=y), color = bg_color, size = point_size) +
    geom_point(data = subset(panel_b_df, display_group != "Other"), 
               aes(x=x, y=y, color = display_group), size = point_size) +
    scale_color_manual(values = panel_b_colors, name = "Fibroblast State") +
    labs(title = "B. Fibroblast State (Annotation)") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3)))

# --- PANEL C: Anomalies (Red/Green - The only place these are used) ---
panel_c_df <- plot_df
panel_c_df$display_group <- "Other"
is_fibro_lineage <- panel_c_df$broad_type %in% c("Fibroblast", "Inflamed Fibroblast")
panel_c_df$display_group[is_fibro_lineage & panel_c_df$anomaly_status == "Typical"] <- "Typical"
panel_c_df$display_group[is_fibro_lineage & panel_c_df$anomaly_status == "Anomalous"] <- "Anomalous"

panel_c_colors <- c(
    "Typical" = "#33A02C",   # Green
    "Anomalous" = "#E31A1C", # Red
    "Other" = bg_color
)

panel_c <- ggplot() +
    geom_point(data = subset(panel_c_df, display_group == "Other"), 
               aes(x=x, y=y), color = bg_color, size = point_size) +
    geom_point(data = subset(panel_c_df, display_group != "Other"), 
               aes(x=x, y=y, color = display_group), size = point_size) +
    scale_color_manual(values = panel_c_colors, name = "scDiagnostics Prediction") +
    labs(title = "C. Detected Anomalies") +
    common_theme +
    guides(color = guide_legend(override.aes = list(size=3)))

# --- PANEL D: ECM Score ---
panel_d <- ggplot() +
    geom_point(data = subset(plot_df, !broad_type %in% c("Fibroblast", "Inflamed Fibroblast")), 
               aes(x=x, y=y), color = bg_color, size = point_size) +
    geom_point(data = subset(plot_df, broad_type %in% c("Fibroblast", "Inflamed Fibroblast")), 
               aes(x=x, y=y, color = ecm_score), size = point_size) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "ECM Score") +
    labs(title = "D. Loss of Homeostasis Signature") +
    common_theme

# -----------------------------------------------------------
# 6. Assemble & Save
# -----------------------------------------------------------
cat("\n--- Assembling Final Figure ---\n")
final_figure <- (panel_a + panel_b) / (panel_c + panel_d) + 
    plot_annotation(
        title = "Spatial Validation of Anomalous Fibroblasts (MERFISH)",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

print(final_figure)

dir.create("figures/merfish", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/merfish/spatial_plots.png", 
       plot = final_figure, width = 14, height = 10, units = "in", dpi = 600)

cat("\n✓ Saved. Check console for Performance Metrics.\n")