# --------------------------------------------------------------------
# AUDIT: scDiagnostics QC of Azimuth Labels (High Specificity Mode)
# + SPATIAL VISUALIZATION + METRICS
# --------------------------------------------------------------------

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scDiagnostics) 

# --- 1. Load Data ---
cat("--- Loading Data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l1 <- as.character(dss9_data$azimuth_celltype_l1)

# --- 2. Align Annotations (LUMPING STRATEGY) ---
# Create temporary columns for the algorithm to use.
# Goal: Compare generic "Fibroblast" Query vs generic "Fibroblast" Reference.

healthy_data$qc_label <- healthy_data$tier2
dss9_data$qc_label    <- dss9_data$azimuth_celltype_l1

# A. Reference Lumping (Healthy)
# We define "Healthy Reference" as any Fibro except Fibro 4 (Inflamed).
# We lump Fibro 1, 2, 3, 5, etc. -> "Fibroblast"
is_healthy_fibro_ref <- grepl("^Fibro", healthy_data$tier2) 
healthy_data$qc_label[is_healthy_fibro_ref] <- "Fibroblast"

# B. Query Lumping (Azimuth)
# If Azimuth calls it "Fibroblast" (or any subtype), we map it to "Fibroblast"
dss9_data$qc_label[grepl("Fibro", dss9_data$azimuth_celltype_l1, ignore.case=TRUE)] <- "Fibroblast"

# Check overlap
common_types <- intersect(unique(dss9_data$qc_label), unique(healthy_data$qc_label))

cat("\n------------------------------------------------------------\n")
cat("Testing Specificity/Sensitivity using Unified Labels.\n")
cat("Reference 'Fibro 1,2,3...' lumped to: 'Fibroblast'\n")
cat("Query Azimuth 'Fibro...' lumped to:    'Fibroblast'\n")
cat("Overlapping Types for Analysis:", paste(common_types, collapse=", "), "\n")
cat("------------------------------------------------------------\n")

if (!"Fibroblast" %in% common_types) {
    stop("Error: 'Fibroblast' category not found in overlap. Check the regex logic.")
}

# --- 3. Run Anomaly Detection (Using QC_LABEL) ---
cat("--- Running Anomaly Detection ---\n")

anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "qc_label",    # <--- Uses the Lumped Label
    ref_cell_type_col = "qc_label",      # <--- Uses the Lumped Label
    cell_types = common_types,                   
    pc_subset = 1:3, 
    anomaly_threshold = 0.5,
    max_cells_query = NULL,
    max_cells_ref = NULL
)

# Extract Anomalous Barcodes
anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# --- 4. The Audit: Focus on Ground Truth Fibroblasts ---
# We still use the original Tier 2 to know what they *actually* are
audit_df <- colData(dss9_data) %>% as.data.frame() %>%
    filter(grepl("^Fibro|^IAF", tier2)) %>%
    mutate(
        Is_Anomalous = ifelse(rownames(.) %in% anomalous_barcodes, TRUE, FALSE),
        Azimuth_Label = azimuth_celltype_l1
    )

# --- 5. Generate Statistics (Grouped by Original Azimuth Label) ---
stats_df <- audit_df %>%
    group_by(Azimuth_Label) %>%
    summarise(
        Total_Real_Fibros = n(), 
        Count_Flagged_Anomalous = sum(Is_Anomalous),
        Percent_Flagged = (sum(Is_Anomalous) / n()) * 100
    ) %>%
    filter(Total_Real_Fibros > 5) %>%
    arrange(desc(Percent_Flagged))

# --- 6. Print Report ---
cat("\n================================================================================\n")
cat(" QC AUDIT REPORT (Lumped Fibroblast Model)\n")
cat("================================================================================\n")
print(as.data.frame(stats_df))

# --- 7. Bar Plot Visualization ---
stats_df$Classification_Type <- ifelse(grepl("Fibro", stats_df$Azimuth_Label, ignore.case = TRUE), 
                                       "Correct Lineage", "Misclassification")

p1 <- ggplot(stats_df, aes(x = reorder(Azimuth_Label, Percent_Flagged), y = Percent_Flagged, fill=Classification_Type)) +
    geom_bar(stat = "identity", color="black", width=0.7) +
    scale_fill_manual(values=c("Correct Lineage" = "grey70", "Misclassification" = "#D9230F")) +
    coord_flip() +
    theme_bw() +
    labs(
        title = "scDiagnostics: Azimuth QC (Unified Fibroblast Model)",
        subtitle = "Anomaly detection performed comparing generic 'Fibroblast' labels.\nMisclassifications (Red) should be high. Healthy Fibros (Grey) should be low.",
        y = "% Flagged as Anomalous",
        x = "Original Azimuth Label",
        fill = "Azimuth Accuracy"
    ) +
    geom_text(aes(label = paste0(round(Percent_Flagged, 0), "% (n=", Total_Real_Fibros, ")")), 
              hjust = -0.1, size=3) +
    ylim(0, 115)

print(p1)
ggsave("figures/merfish/qc_audit_azimuth_unified_bars.png", plot=p1, width=9, height=8)

# --- 8. SPATIAL VISUALIZATION (Confusion Matrix Map) ---
cat("\n--- Generating Spatial Map of Sensitivity/Specificity ---\n")

# A. Prepare Data
spatial_df <- colData(dss9_data) %>% as.data.frame()

# B. Coordinate Finder
if("center_x" %in% colnames(spatial_df) & "center_y" %in% colnames(spatial_df)) {
    spatial_df$plot_x <- spatial_df$center_x
    spatial_df$plot_y <- spatial_df$center_y
} else if("x" %in% colnames(spatial_df) & "y" %in% colnames(spatial_df)) {
    spatial_df$plot_x <- spatial_df$x
    spatial_df$plot_y <- spatial_df$y
} else {
    coords <- spatialCoords(dss9_data)
    spatial_df$plot_x <- coords[,1]
    spatial_df$plot_y <- coords[,2]
}

# C. Define Categories based on Ground Truth & Prediction
# Filter to only Fibroblast Lineage (Background cells get a separate tag)
spatial_df$Map_Status <- "Background (Non-Fibro)"
is_fibro_lineage <- grepl("^Fibro|^IAF", spatial_df$tier2)

# 1. Prediction: Did scDiagnostics flag it?
is_anomalous <- rownames(spatial_df) %in% anomalous_barcodes

# 2. Truth: Is it Inflamed? (IAF or Fibro 4)
is_inflamed <- grepl("^IAF", spatial_df$tier2)

# 3. Assign Confusion Matrix Categories
# Only apply this logic to the Ground Truth Fibroblasts
spatial_df$Map_Status[is_fibro_lineage] <- case_when(
    # TRUE POSITIVE: Inflamed AND Flagged Anomalous
    is_inflamed[is_fibro_lineage] & is_anomalous[is_fibro_lineage] ~ "TP: Inflamed Detected (Success)",
    
    # FALSE NEGATIVE: Inflamed BUT Called Typical
    is_inflamed[is_fibro_lineage] & !is_anomalous[is_fibro_lineage] ~ "FN: Inflamed Missed",
    
    # TRUE NEGATIVE: Healthy AND Called Typical
    !is_inflamed[is_fibro_lineage] & !is_anomalous[is_fibro_lineage] ~ "TN: Healthy Typical",
    
    # FALSE POSITIVE: Healthy BUT Flagged Anomalous
    !is_inflamed[is_fibro_lineage] & is_anomalous[is_fibro_lineage] ~ "FP: Healthy Flagged (Noise)"
)

# D. Set Plotting Order (Prioritize TP/FN/FP over TN)
spatial_df <- spatial_df %>%
    mutate(order_val = case_when(
        Map_Status == "Background (Non-Fibro)" ~ 1,
        Map_Status == "TN: Healthy Typical" ~ 2,
        Map_Status == "FN: Inflamed Missed" ~ 3,
        Map_Status == "FP: Healthy Flagged (Noise)" ~ 4,
        Map_Status == "TP: Inflamed Detected (Success)" ~ 5
    )) %>%
    arrange(order_val)

# E. Define Colors
map_colors <- c(
    "Background (Non-Fibro)"          = "grey95",  
    "TN: Healthy Typical"             = "grey70",  # Matches previous 'Correct'
    "TP: Inflamed Detected (Success)" = "#D9230F", # RED (Sensitivity)
    "FN: Inflamed Missed"             = "black",   # BLACK (Missed sensitivity)
    "FP: Healthy Flagged (Noise)"     = "#E69F00"  # ORANGE (Specificity hit)
)

# F. Generate Plot
p2 <- ggplot(spatial_df, aes(x = plot_x, y = plot_y, color = Map_Status)) +
    geom_point(size = 0.8, stroke = 0) +
    scale_color_manual(values = map_colors) +
    theme_void() +
    coord_fixed() + 
    guides(color = guide_legend(override.aes = list(size=3))) +
    labs(
        title = "Spatial Map: Detection of Inflamed States",
        subtitle = "Performance of anomaly detection on Unified Fibroblast Labels",
        color = "Performance Class"
    )

print(p2)
ggsave("figures/merfish/qc_audit_spatial_sensitivity.png", plot=p2, width=10, height=8, dpi=600)

# --- 9. SENSITIVITY & SPECIFICITY METRICS ---
cat("\n================================================================================\n")
cat(" PERFORMANCE METRICS: Inflamed vs Healthy Fibro Detection\n")
cat(" Method: Unified 'Fibroblast' Model (Healthy Fibros lumped)\n")
cat("================================================================================\n")

# 1. Define Population & Truth
perf_df <- colData(dss9_data) %>% as.data.frame() %>%
    filter(grepl("^Fibro|^IAF", tier2)) %>%
    mutate(
        Predicted_Anomalous = rownames(.) %in% anomalous_barcodes,
        # Truth: Inflamed = IAF subtypes OR Fibro 4
        Is_Inflamed_Truth = grepl("IAF", tier2) 
    )

# 2. Confusion Matrix
TP <- sum(perf_df$Predicted_Anomalous == TRUE & perf_df$Is_Inflamed_Truth == TRUE)
FN <- sum(perf_df$Predicted_Anomalous == FALSE & perf_df$Is_Inflamed_Truth == TRUE)
TN <- sum(perf_df$Predicted_Anomalous == FALSE & perf_df$Is_Inflamed_Truth == FALSE)
FP <- sum(perf_df$Predicted_Anomalous == TRUE & perf_df$Is_Inflamed_Truth == FALSE)

# 3. Metrics
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
accuracy    <- (TP + TN) / (TP + FN + TN + FP)

# 4. Report
cat(sprintf("\nConfusion Matrix (n=%d Total Fibros):\n", nrow(perf_df)))
cat(sprintf("                 | True Inflamed | True Healthy\n"))
cat(sprintf("  Pred Anomalous | %-13d | %-12d (TP | FP)\n", TP, FP))
cat(sprintf("  Pred Typical   | %-13d | %-12d (FN | TN)\n\n", FN, TN))

cat(sprintf("  Sensitivity (Detect Inflamed/IAF):  %.2f%%\n", sensitivity * 100))
cat(sprintf("  Specificity (Keep Healthy Typical): %.2f%%\n", specificity * 100))
cat(sprintf("  Overall Accuracy:                   %.2f%%\n", accuracy * 100))
cat("================================================================================\n")