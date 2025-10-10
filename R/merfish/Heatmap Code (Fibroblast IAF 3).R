# -----------------------------------------------
# MERFISH Mouse Colon IBD - IAF3 Problem vs. Solution Figures 
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ComplexHeatmap)
library(circlize)
library(stringr)

# -----------------------------------------------

# Load and prepare data
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l2 <- as.character(dss9_data$azimuth_celltype_l2) 

# ... [The Anomaly Detection part is correct and remains the same] ...
cat("--- Running anomaly detection on naive 'Fibro 2' and 'Fibro 6' predictions ---\n")
anomaly_output <- detectAnomaly(query_data = dss9_data, reference_data = healthy_data, 
    query_cell_type_col = "azimuth_celltype_l2", ref_cell_type_col = "tier2", 
    cell_types = c("Fibro 2", "Fibro 6"), pc_subset = 1:5,
    anomaly_threshold = 0.5)
anomalous_barcodes <- gsub("Query_", "", c(names(anomaly_output$`Fibro 2`$query_anomaly[anomaly_output$`Fibro 2`$query_anomaly == TRUE]), names(anomaly_output$`Fibro 6`$query_anomaly[anomaly_output$`Fibro 6`$query_anomaly == TRUE])))
dss9_data$anomaly_status_new <- "Typical"
dss9_data$anomaly_status_new[colnames(dss9_data) %in% anomalous_barcodes] <- "Anomalous"
cat("Anomaly detection complete.\n\n")

# --- Define Shared Variables ---
genes_fibro6 <- c("Wnt2b", "Grem1", "Pdgfra")
genes_fibro2 <- c("Wnt5a", "Adamdec1", "Bmp5")
genes_iaf3 <- c("Igfbp5", "Col18a1", "Mmp2", "Cxcl13")
genes_of_interest <- unique(c(genes_fibro6, genes_fibro2, genes_iaf3))


####################################################################
# FIGURE A: The "Problem" - Mislabeled Groups Obscure Biology
####################################################################
cat("--- Generating Figure A: The 'Problem' Heatmap ---\n")

# --- Define groups for Figure A ---
groups_a <- list()
groups_a$`Healthy Fibro 2` <- which(healthy_data$tier2 == "Fibro 2")
groups_a$`Healthy Fibro 6` <- which(healthy_data$tier2 == "Fibro 6")
groups_a$`DSS9 "Fibro 2" (Azimuth)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 2")
groups_a$`DSS9 "Fibro 6" (Azimuth)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 6")
groups_a$`DSS9 IAF3 (Ground Truth)` <- which(dss9_data$tier2 == "IAF 3")

# --- Calculate and Scale Matrix for Figure A ---
# ... [This logic is correct and remains the same] ...
genes_present_a <- genes_of_interest[genes_of_interest %in% rownames(healthy_data)]
gene_indices_healthy_a <- match(genes_present_a, rownames(healthy_data))
gene_indices_dss9_a <- match(genes_present_a, rownames(dss9_data))
pseudo_bulk_a <- sapply(groups_a, function(indices) {
    if (length(indices) == 0) return(rep(NA, length(genes_present_a)))
    sce_obj <- if (any(sapply(names(indices), function(n) n %in% colnames(healthy_data)))) healthy_data else dss9_data
    gene_idx <- if (identical(sce_obj, healthy_data)) gene_indices_healthy_a else gene_indices_dss9_a
    rowMeans(as.matrix(logcounts(sce_obj)[gene_idx, indices]))
})
rownames(pseudo_bulk_a) <- genes_present_a
final_matrix_a <- pseudo_bulk_a[, colSums(is.na(pseudo_bulk_a)) == 0]
scaled_matrix_a <- t(scale(t(final_matrix_a)))

# --- Create Annotations for Figure A (CORRECTED) ---
ann_df_a <- data.frame(row.names = colnames(scaled_matrix_a))
ann_df_a$Source <- ifelse(grepl("Healthy", colnames(scaled_matrix_a)), "Reference", "Query")

# *** KEY FIX IS HERE: Use str_extract for precise matching ***
ann_df_a$`Cell Type` <- str_extract(colnames(scaled_matrix_a), "Fibro 2|Fibro 6|IAF3")

ann_colors_a <- list(
    Source = c(Reference = "navy", Query = "darkred"),
    `Cell Type` = c(`Fibro 2` = "steelblue", `Fibro 6` = "skyblue", `IAF3` = "firebrick")
)
top_ann_a <- HeatmapAnnotation(df = ann_df_a, col = ann_colors_a, annotation_name_side = "left")

# --- Draw Figure A ---
ht_a <- Heatmap(scaled_matrix_a, name = "Scaled Expression", col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
                top_annotation = top_ann_a, cluster_columns = FALSE, column_order = names(groups_a)[names(groups_a) %in% colnames(scaled_matrix_a)],
                column_title = "Figure A: Biology is Obscured by Mislabeled Cells",
                column_names_rot = 45, column_names_gp = gpar(fontsize=10))
draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right")


####################################################################
# FIGURE B: The "Solution" - Anomaly Detection Recovers Biology
####################################################################
cat("\n--- Generating Figure B: The 'Solution' Heatmap ---\n")

# --- Define groups for Figure B ---
groups_b <- list()
groups_b$`Healthy Fibro 2` <- which(healthy_data$tier2 == "Fibro 2")
groups_b$`Healthy Fibro 6` <- which(healthy_data$tier2 == "Fibro 6")
groups_b$`DSS9 "Fibro 2" (Typical)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 2" & dss9_data$anomaly_status_new == "Typical")
groups_b$`DSS9 "Fibro 2" (Anomalous)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 2" & dss9_data$anomaly_status_new == "Anomalous")
groups_b$`DSS9 "Fibro 6" (Typical)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 6" & dss9_data$anomaly_status_new == "Typical")
groups_b$`DSS9 "Fibro 6" (Anomalous)` <- which(dss9_data$azimuth_celltype_l2 == "Fibro 6" & dss9_data$anomaly_status_new == "Anomalous")
groups_b$`DSS9 IAF3 (Ground Truth)` <- which(dss9_data$tier2 == "IAF 3")

# --- Calculate and Scale Matrix for Figure B ---
# ... [This logic is correct and remains the same] ...
genes_present_b <- genes_of_interest[genes_of_interest %in% rownames(healthy_data)]
gene_indices_healthy_b <- match(genes_present_b, rownames(healthy_data))
gene_indices_dss9_b <- match(genes_present_b, rownames(dss9_data))
pseudo_bulk_b <- sapply(groups_b, function(indices) {
    if (length(indices) == 0) return(rep(NA, length(genes_present_b)))
    sce_obj <- if (any(sapply(names(indices), function(n) n %in% colnames(healthy_data)))) healthy_data else dss9_data
    gene_idx <- if (identical(sce_obj, healthy_data)) gene_indices_healthy_b else gene_indices_dss9_b
    rowMeans(as.matrix(logcounts(sce_obj)[gene_idx, indices]))
})
rownames(pseudo_bulk_b) <- genes_present_b
final_matrix_b <- pseudo_bulk_b[, colSums(is.na(pseudo_bulk_b)) == 0]
scaled_matrix_b <- t(scale(t(final_matrix_b)))


# --- Create Annotations for Figure B (CORRECTED) ---
ann_df_b <- data.frame(row.names = colnames(scaled_matrix_b))
ann_df_b$Source <- ifelse(grepl("Healthy", colnames(scaled_matrix_b)), "Reference", "Query")

# *** KEY FIX IS HERE: Use str_extract for precise matching ***
ann_df_b$`Cell Type` <- str_extract(colnames(scaled_matrix_b), "Fibro 2|Fibro 6|IAF3")

ann_df_b$`Anomaly Status` <- "N/A"
ann_df_b$`Anomaly Status`[grepl("Typical", colnames(scaled_matrix_b))] <- "Typical"
ann_df_b$`Anomaly Status`[grepl("Anomalous", colnames(scaled_matrix_b))] <- "Anomalous"

ann_colors_b <- list(
    Source = c(Reference = "navy", Query = "darkred"),
    `Cell Type` = c(`Fibro 2` = "steelblue", `Fibro 6` = "skyblue", `IAF3` = "firebrick"),
    `Anomaly Status` = c(Typical = "grey80", Anomalous = "gold", "N/A" = "white")
)
top_ann_b <- HeatmapAnnotation(df = ann_df_b, col = ann_colors_b, annotation_name_side = "left")

# --- Draw Figure B ---
ht_b <- Heatmap(scaled_matrix_b, name = "Scaled Expression", col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
                top_annotation = top_ann_b, cluster_columns = FALSE, column_order = names(groups_b)[names(groups_b) %in% colnames(scaled_matrix_b)],
                column_title = "Figure B: Anomaly Detection Recovers True Biological States",
                column_names_rot = 45, column_names_gp = gpar(fontsize=10))
draw(ht_b, heatmap_legend_side = "right", annotation_legend_side = "right")