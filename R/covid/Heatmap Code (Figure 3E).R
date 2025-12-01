# --------------------------------------
# COVID-19 - Heatmaps Code (Figure 3E)
# --------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(scDiagnostics) 

# --------------------------------------

# Load Data
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# ____________________________
# 1. Detect Anomalies
# ____________________________
# Run detection for both cell types. This is essential.
anomaly_output <- detectAnomaly(
    query_data = covid_data,
    reference_data = normal_data, 
    query_cell_type_col = "singler_annotations", 
    ref_cell_type_col = "cell_type", 
    cell_types = c("CD14-positive monocyte", "mature NK T cell"), 
    pc_subset = 1:10
)

# _________________
# 2. Add Anomaly Status to Data
# _________________
anomalous_mono_barcodes <- names(anomaly_output$`CD14-positive monocyte`$query_anomaly[anomaly_output$`CD14-positive monocyte`$query_anomaly == TRUE])
anomalous_nkt_barcodes <- names(anomaly_output$`mature NK T cell`$query_anomaly[anomaly_output$`mature NK T cell`$query_anomaly == TRUE])
all_anomalous_barcodes <- unique(c(anomalous_mono_barcodes, anomalous_nkt_barcodes))
covid_data$anomaly_status <- "Typical"
covid_data$anomaly_status[colnames(covid_data) %in% all_anomalous_barcodes] <- "Anomalous"

# __________________________________________________
# 3. Build and Scale Data for EACH Heatmap SEPARATELY
# __________________________________________________

# Define global variables
severity_levels <- c("Asymptomatic", "Mild", "Moderate", "Severe", "Critical")
genes_mono <- c("CD14", "LYZ", "VCAN")
genes_nkt <- c("NKG7", "KLRD1", "CD3D", "GZMB")
genes_inflam <- c("S100A8", "S100A9", "S100A12", "ISG15", "MX1", "IFITM3")

# --- Part A: Build the Monocyte Matrix ---
mono_genes_for_plot <- unique(c(genes_mono, genes_inflam))
mono_genes_present <- mono_genes_for_plot[mono_genes_for_plot %in% rowData(normal_data)$gene_symbol]
mono_gene_indices <- match(mono_genes_present, rowData(normal_data)$gene_symbol)

mono_groups <- list()
mono_groups$`Monocyte (Healthy)` <- which(normal_data$cell_type == "CD14-positive monocyte")
for (sev in severity_levels) {
    mono_groups[[paste0("Monocyte (", sev, ", Typical)")]] <- which(covid_data$cell_type == "CD14-positive monocyte" & covid_data$Status_on_day_collection_summary == sev & covid_data$anomaly_status == "Typical")
    mono_groups[[paste0("Monocyte (", sev, ", Anomalous)")]] <- which(covid_data$cell_type == "CD14-positive monocyte" & covid_data$Status_on_day_collection_summary == sev & covid_data$anomaly_status == "Anomalous")
}
mono_pseudo_bulk <- sapply(mono_groups, function(indices) {
    if (length(indices) == 0) return(rep(NA, length(mono_genes_present)))
    is_normal <- any(sapply(names(indices), function(n) n %in% colnames(normal_data)))
    sce_object <- if(is_normal) normal_data else covid_data
    rowMeans(as.matrix(logcounts(sce_object)[mono_gene_indices, indices, drop = FALSE]))
})
rownames(mono_pseudo_bulk) <- mono_genes_present
mono_final_matrix <- mono_pseudo_bulk[, colSums(is.na(mono_pseudo_bulk)) == 0, drop = FALSE]
mono_plot_matrix <- t(scale(t(mono_final_matrix)))

# --- Part B: Build the NK T Cell Matrix (CORRECTED) ---
nkt_genes_for_plot <- unique(c(genes_nkt, genes_inflam))
# *** KEY FIX: Use the same inflammatory genes that actually made it into the mono matrix ***
mono_inflam_genes_present <- rownames(mono_plot_matrix)[rownames(mono_plot_matrix) %in% genes_inflam]
nkt_genes_for_plot <- unique(c(genes_nkt, mono_inflam_genes_present))  # Use SAME inflam genes

nkt_genes_present <- nkt_genes_for_plot[nkt_genes_for_plot %in% rowData(normal_data)$gene_symbol]
nkt_gene_indices <- match(nkt_genes_present, rowData(normal_data)$gene_symbol)

nkt_groups <- list()
for (sev in severity_levels) {
    nkt_groups[[paste0("NK T (", sev, ", Correct)")]] <- which(covid_data$cell_type == "mature NK T cell" & covid_data$singler_annotations == "mature NK T cell" & covid_data$Status_on_day_collection_summary == sev)
    nkt_groups[[paste0("NK T (", sev, ", Rescued)")]] <- which(covid_data$cell_type == "mature NK T cell" & covid_data$singler_annotations == "CD14-positive monocyte" & covid_data$anomaly_status == "Anomalous" & covid_data$Status_on_day_collection_summary == sev)
    nkt_groups[[paste0("NK T (", sev, ", Missed)")]] <- which(covid_data$cell_type == "mature NK T cell" & covid_data$singler_annotations == "CD14-positive monocyte" & covid_data$anomaly_status == "Typical" & covid_data$Status_on_day_collection_summary == sev)
}
nkt_pseudo_bulk <- sapply(nkt_groups, function(indices) {
    if (length(indices) == 0) return(rep(NA, length(nkt_genes_present)))
    rowMeans(as.matrix(logcounts(covid_data)[nkt_gene_indices, indices, drop = FALSE]))
})
rownames(nkt_pseudo_bulk) <- nkt_genes_present
nkt_final_matrix <- nkt_pseudo_bulk[, colSums(is.na(nkt_pseudo_bulk)) == 0, drop = FALSE]
nkt_plot_matrix <- t(scale(t(nkt_final_matrix)))

# *** ALSO UPDATE: Create the NK T row split using the same genes ***
nkt_row_split <- factor(ifelse(rownames(nkt_plot_matrix) %in% mono_inflam_genes_present, "Inflammatory Markers", "NK T Markers"), levels = c("NK T Markers", "Inflammatory Markers"))

# _________________________________________________________________
# 4. Prepare Annotations and Plotting Variables
# _________________________________________________________________

# Create a single, master annotation object for ALL columns that will appear
all_cols <- c(colnames(mono_plot_matrix), colnames(nkt_plot_matrix))
col_data_for_annotation <- data.frame(row.names = all_cols)
col_data_for_annotation$`Ground Truth` <- ifelse(grepl("Monocyte", all_cols), "Monocyte", "NK T")
col_data_for_annotation$`Package Output` <- "Typical"
col_data_for_annotation$`Package Output`[grepl("Anomalous|Rescued", all_cols)] <- "Anomalous"
col_data_for_annotation$`Package Output`[grepl("Missed", all_cols)] <- "Missed (Anomalous)"
col_data_for_annotation$`Disease Status` <- "Healthy"
for (sev in severity_levels) { col_data_for_annotation$`Disease Status`[grepl(sev, all_cols)] <- sev }

ann_colors <- list( `Ground Truth` = c("Monocyte" = "#1F78B4", "NK T" = "#33A02C"), `Package Output` = c("Typical" = "grey80", "Anomalous" = "#E31A1C", "Missed (Anomalous)" = "#FB9A99"), `Disease Status` = c("Healthy" = "white", "Asymptomatic" = "#FFFFB2", "Mild" = "#FED976", "Moderate" = "#FEB24C", "Severe" = "#FD8D3C", "Critical" = "#E31A1C") )
top_annotation <- HeatmapAnnotation(df = col_data_for_annotation, col = ann_colors, annotation_name_side = "left", annotation_name_gp = gpar(fontsize=10))

# Create robust row splits for each panel
mono_row_split <- factor(ifelse(rownames(mono_plot_matrix) %in% genes_inflam, "Inflammatory Markers", "Monocyte Markers"), levels = c("Monocyte Markers", "Inflammatory Markers"))
nkt_row_split <- factor(ifelse(rownames(nkt_plot_matrix) %in% genes_inflam, "Inflammatory Markers", "NK T Markers"), levels = c("NK T Markers", "Inflammatory Markers"))

# _________________________________________________________________
# 5. Define and Draw the Composite Heatmap (CORRECTED)
# _________________________________________________________________

# Create separate annotation objects for each heatmap
mono_col_data <- col_data_for_annotation[colnames(mono_plot_matrix), ]
mono_annotation <- HeatmapAnnotation(df = mono_col_data, col = ann_colors, 
                                     annotation_name_side = "left", 
                                     annotation_name_gp = gpar(fontsize=10))

nkt_col_data <- col_data_for_annotation[colnames(nkt_plot_matrix), ]
nkt_annotation <- HeatmapAnnotation(df = nkt_col_data, col = ann_colors, 
                                    annotation_name_side = "left", 
                                    annotation_name_gp = gpar(fontsize=10))

# --- LEFT heatmap: Monocytes ---
ht_mono <- Heatmap(
    mono_plot_matrix, name = "Scaled Expression", col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
    top_annotation = mono_annotation,  # Use the new annotation object
    cluster_columns = FALSE,
    column_names_rot = 90, column_names_gp = gpar(fontsize = 8),
    row_split = mono_row_split, cluster_rows = TRUE, 
    row_gap = unit(3, "mm"),  # We can add this back now
    row_names_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10, fontface = "italic"), row_title_rot = 0,
    show_heatmap_legend = FALSE
)

# --- RIGHT heatmap: NK T Cells ---
ht_nkt <- Heatmap(
    nkt_plot_matrix, name = "Scaled Expression", col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
    top_annotation = nkt_annotation,  # Use the new annotation object
    cluster_columns = FALSE, column_title = "NK T Cell Deconvolution",
    column_names_rot = 90, column_names_gp = gpar(fontsize = 8),
    row_split = nkt_row_split, cluster_rows = TRUE,
    row_gap = unit(3, "mm"),  # We can add this back now
    show_row_names = FALSE, row_title = NULL, show_heatmap_legend = FALSE
)
