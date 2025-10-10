# -----------------------------------------------
# MERFISH Mouse Colon IBD - Fibro 7/4 Figure (Final with Annotations)
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics) 
library(ComplexHeatmap)
library(circlize)

# -----------------------------------------------

# Load and prepare data
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l2 <- as.character(dss9_data$azimuth_celltype_l2)

# __________________________________________________________________
# 1. Run Anomaly Detection
# __________________________________________________________________
cat("Running anomaly detection on the naive Azimuth 'Fibro 7' predictions...\n")
anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "azimuth_celltype_l2",
    ref_cell_type_col = "tier2", 
    cell_types = c("Fibro 7"), 
    pc_subset = 1:10
)
anomalous_barcodes <- gsub("Query_", "", names(anomaly_output$`Fibro 7`$query_anomaly[anomaly_output$`Fibro 7`$query_anomaly == TRUE]))
dss9_data$anomaly_status_fibro7 <- "Typical"
dss9_data$anomaly_status_fibro7[colnames(dss9_data) %in% anomalous_barcodes] <- "Anomalous"

# ______________________________________________________
# 2. Prepare Data and Annotations for the Heatmap
# ______________________________________________________

# --- Gene list and cell groups (same as before) ---
genes_fibro7_id <- c("Cd34", "Pi16", "Ackr4")
genes_myofibro <- c("Tagln", "Dpt", "Mmp3")
genes_inflam_sig <- c("Il1b", "Il11", "Cxcl12", "Csf3")
genes_of_interest <- unique(c(genes_fibro7_id, genes_myofibro, genes_inflam_sig))

group_indices <- list()
all_cells_list <- list() # Store cell indices for validation bar plot

# Panel A Groups
group_indices$`Healthy Fibro 7 (A)` <- which(healthy_data$tier2 == "Fibro 7")
group_indices$`DSS9 "Fibro 7" Naive (A)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 7")
group_indices$`Ground Truth Fibro 4 (A)` <- which(dss9_data$tier2 == "Fibro 4")

# Panel B Groups
group_indices$`Healthy Fibro 7 (B)` <- which(healthy_data$tier2 == "Fibro 7")
group_indices$`DSS9 "Fibro 7" Typical (B)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 7" & dss9_data$anomaly_status_fibro7 == "Typical")
group_indices$`DSS9 "Fibro 7" Anomalous (B)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 7" & dss9_data$anomaly_status_fibro7 == "Anomalous")
group_indices$`Ground Truth Fibro 4 (B)` <- which(dss9_data$tier2 == "Fibro 4")

# --- Calculate and Scale the Matrix (same as before) ---
# ... [This logic is robust and unchanged] ...
genes_present <- genes_of_interest[genes_of_interest %in% rownames(healthy_data)]
gene_indices_healthy <- match(genes_present, rownames(healthy_data))
gene_indices_dss9 <- match(genes_present, rownames(dss9_data))
pseudo_bulk_list <- list()
for (group_name in names(group_indices)) {
    indices <- group_indices[[group_name]]
    all_cells_list[[group_name]] <- if (grepl("Healthy", group_name)) healthy_data[, indices] else dss9_data[, indices] # Store the sce subset
    if (length(indices) == 0) { pseudo_bulk_list[[group_name]] <- rep(NA, length(genes_present)); next }
    sce_obj <- if (grepl("Healthy", group_name)) healthy_data else dss9_data
    gene_idx <- if (grepl("Healthy", group_name)) gene_indices_healthy else gene_indices_dss9
    pseudo_bulk_list[[group_name]] <- rowMeans(as.matrix(logcounts(sce_obj)[gene_idx, indices]))
}
pseudo_bulk_matrix <- do.call(cbind, pseudo_bulk_list)
final_matrix <- pseudo_bulk_matrix[, colSums(is.na(pseudo_bulk_matrix)) == 0]
scaled_matrix <- t(scale(t(final_matrix)))
rownames(scaled_matrix) <- genes_present[genes_present %in% rownames(final_matrix)]


# --- *** NEW: Create the Multi-Layered Annotation Object *** ---

# Clean up column names for the plot
plot_column_names <- gsub(" \\(A\\)| \\(B\\)", "", colnames(scaled_matrix))
names(plot_column_names) <- colnames(scaled_matrix)

# Annotation 1: Simple Group Labels (for both panels)
group_colors <- c(
    `Healthy Fibro 7` = "steelblue", `DSS9 "Fibro 7" Naive` = "grey50", `Ground Truth Fibro 4` = "firebrick",
    `DSS9 "Fibro 7" Typical` = "lightblue", `DSS9 "Fibro 7" Anomalous` = "darkorange"
)
group_annotation_vector <- plot_column_names

# Annotation 2: Validation Bar Plot (for Panel B only)
composition_list <- lapply(names(all_cells_list), function(group_name) {
    sce_subset <- all_cells_list[[group_name]]
    # Only calculate for columns that exist in our final matrix
    if (!group_name %in% colnames(scaled_matrix)) return(NA)
    is_fibro4 <- sce_subset$tier2 == "Fibro 4"
    if (length(is_fibro4) == 0) return(0)
    round(100 * sum(is_fibro4) / length(is_fibro4))
})
names(composition_list) <- names(all_cells_list)
validation_barplot_vector <- unlist(composition_list)
validation_barplot_vector <- validation_barplot_vector[colnames(scaled_matrix)] # Ensure order and presence
# For Panel A, we want this to be empty. We set their values to NA.
validation_barplot_vector[grepl("\\(A\\)", names(validation_barplot_vector))] <- NA


# Create the final HeatmapAnnotation object
top_annotation <- HeatmapAnnotation(
    # Layer 1: Colored rectangles for group identity
    Group = anno_simple(group_annotation_vector, col = group_colors),
    
    # Layer 2: Bar plot for validation
    `% Fibro 4` = anno_barplot(validation_barplot_vector, 
                               gp = gpar(fill = "black"),
                               ylim = c(0, 100),
                               axis_param = list(gp = gpar(fontsize = 8))),
    
    # Show the bar plot only for the right-hand panel (Panel B)
    which = "column",
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize=10)
)


# _________________________________________________________________
# 3. Create and Draw the Final, Single Heatmap
# _________________________________________________________________

# --- Define row and column splitting factors ---
row_splitting_factor <- factor(
    ifelse(rownames(scaled_matrix) %in% genes_inflam_sig, "Pro-inflammatory Signaling", 
           ifelse(rownames(scaled_matrix) %in% genes_myofibro, "Myofibroblast / ECM", "Homeostatic Fibro 7 Identity")),
    levels = c("Homeostatic Fibro 7 Identity", "Myofibroblast / ECM", "Pro-inflammatory Signaling")
)
column_splitting_factor <- factor(
  ifelse(grepl("\\(A\\)", colnames(scaled_matrix)), "A) Before Anomaly Detection", "B) After Anomaly Detection"),
  levels = c("A) Before Anomaly Detection", "B) After Anomaly Detection")
)

# --- Define final column order ---
final_column_order <- c( 'Healthy Fibro 7 (A)', 'DSS9 "Fibro 7" Naive (A)', 'Ground Truth Fibro 4 (A)',
                         'Healthy Fibro 7 (B)', 'DSS9 "Fibro 7" Typical (B)', 'DSS9 "Fibro 7" Anomalous (B)', 'Ground Truth Fibro 4 (B)' )
final_column_order <- final_column_order[final_column_order %in% colnames(scaled_matrix)]

# --- Draw the single heatmap ---
ht_final <- Heatmap(
    scaled_matrix, 
    name = "Scaled Expression",
    col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
    
    top_annotation = top_annotation,
    
    column_split = column_splitting_factor,
    cluster_column_slices = FALSE,
    column_gap = unit(5, "mm"),
    
    cluster_columns = FALSE,
    column_order = final_column_order,
    column_labels = plot_column_names[final_column_order],
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize=12, fontface="bold"),

    row_split = row_splitting_factor,
    cluster_rows = TRUE,
    row_gap = unit(3, "mm"),
    row_title_gp = gpar(fontsize = 10, fontface = "italic"),
    row_title_rot = 0,
    row_names_gp = gpar(fontsize = 10)
)

# Draw the plot with legends on the right
draw(
    ht_final, 
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(2, 2, 2, 2), "mm")
)

# --- Save the plot ---
pdf("MERFISH_Fibro7_MultiAnnotation_Figure.pdf", width = 12, height = 8)
draw(ht_final, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()