# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 3 Code
# -----------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(scDiagnostics)
library(GGally)
library(ggridges)
library(viridis)
library(cowplot)

# -----------------------------------------------

# __________
# Load data
# __________

normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

set.seed(0)

# ______________________
# Setup and Parameters
# ______________________

yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

CELL_TYPES_TO_PLOT <- c("CD14 mono")
PC_SUBSET <- 1:3
ASSAY_NAME <- "logcounts"
MAX_CELLS <- 1250 # NULL for all cells

# Find available signature genes
common_genes <- intersect(rownames(covid_data), rownames(normal_data))
signature_genes_avail <- intersect(yoshida_ifn_signature, common_genes)

# PC plot names with variance explained
pc_plot_names <- paste0(
    "PC", PC_SUBSET, " (",
    sprintf("%.1f%%", attributes(reducedDim(normal_data, "PCA"))[["percentVar"]][PC_SUBSET]), ")"
)

# ______________________
# Run Anomaly Detection
# ______________________

# SingleR
anomaly_singler <- detectAnomaly(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "singler_annotations_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    n_tree = 500,
    anomaly_threshold = 0.5,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# Azimuth
anomaly_azimuth <- detectAnomaly(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "azimuth_celltype_l1_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    n_tree = 500,
    anomaly_threshold = 0.5,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# CellTypist
anomaly_celltypist <- detectAnomaly(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "celltypist_predicted_labels_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    n_tree = 500,
    anomaly_threshold = 0.5,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# scArches
anomaly_scarches <- detectAnomaly(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "scvi_prediction_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    n_tree = 500,
    anomaly_threshold = 0.5,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# Extract anomaly logical vectors
singler_anomaly_logical <- anomaly_singler[["CD14 mono"]][["query_anomaly"]]
names(singler_anomaly_logical) <- gsub("Query_", "", names(singler_anomaly_logical))

azimuth_anomaly_logical <- anomaly_azimuth[["CD14 mono"]][["query_anomaly"]]
names(azimuth_anomaly_logical) <- gsub("Query_", "", names(azimuth_anomaly_logical))

celltypist_anomaly_logical <- anomaly_celltypist[["CD14 mono"]][["query_anomaly"]]
names(celltypist_anomaly_logical) <- gsub("Query_", "", names(celltypist_anomaly_logical))

scarches_anomaly_logical <- anomaly_scarches[["CD14 mono"]][["query_anomaly"]]
names(scarches_anomaly_logical) <- gsub("Query_", "", names(scarches_anomaly_logical))

# ______________________
# Define Plot Functions
# ______________________

scatter_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
        geom_point(data = ~subset(., dataset == "Query"), alpha = 0.5, size = 1.5, show.legend = FALSE) +
        geom_point(data = ~subset(., dataset == "Reference"), alpha = 0.7, size = 1.5, show.legend = FALSE) +
        scale_shape_manual(values = c("Reference" = 1, "Query" = 16)) +
        viridis::scale_color_viridis(option = "B") +
        theme_minimal() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

ridge_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = aes(x = !!mapping$x, y = dataset, fill = dataset, color = dataset)) +
        ggridges::geom_density_ridges(alpha = 0.6, scale = 1.5, rel_min_height = 0.01) +
        scale_fill_manual(values = c("Reference" = "#5A9BD8", "Query" = "#B565D8"), guide = "none") +
        scale_color_manual(values = c("Reference" = "#5A9BD8", "Query" = "#B565D8"), guide = "none") +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
}

blank_fn <- function(data, mapping, ...) {
    ggplot() + theme_void() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

# _________________________________
# Function to create PCA plot grid
# _________________________________

create_pca_plot_grid <- function(pca_data, color_var, title, method_name) {
    
    # Create legend plot for this color variable
    legend_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = !!sym(color_var), shape = dataset)) +
        geom_point() +
        scale_shape_manual(name = "Dataset", values = c("Reference" = 1, "Query" = 16)) +
        viridis::scale_color_viridis(option = "B", name = title) +
        theme(
            legend.position = "right", 
            legend.box = "vertical",
            legend.key = element_rect(fill = "white", color = NA)
        )
    
    plot_legend <- GGally::grab_legend(legend_plot)
    
    # Create ggpairs plot
    pca_plot <- GGally::ggpairs(
        pca_data,
        columns = paste0("PC", PC_SUBSET),
        columnLabels = pc_plot_names,
        mapping = aes(color = !!sym(color_var), shape = dataset),
        lower = list(continuous = scatter_fn),
        diag = list(continuous = ridge_fn),
        upper = list(continuous = blank_fn),
        progress = FALSE,
        legend = plot_legend,
        title = paste0(method_name, ": ", title, " (CD14+ Monocytes)")
    ) +
    theme(
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black", size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937")
    )
    
    return(pca_plot)
}

# __________________________________
# Function to Add Annotation Scores
# __________________________________

add_scores_to_pca <- function(pca_data, covid_data, normal_data, score_col_name, query_col_name) {
    
    # Separate reference and query indices
    ref_mask <- pca_data$dataset == "Reference"
    query_mask <- pca_data$dataset == "Query"
    
    # Initialize score vector
    scores <- rep(NA, nrow(pca_data))
    
    # Extract original cell IDs from rownames
    original_ids <- gsub("Reference_|Query_", "", rownames(pca_data))
    
    # Match query cells to covid_data
    query_ids <- original_ids[query_mask]
    query_cell_matches <- match(query_ids, colnames(covid_data))
    scores[query_mask] <- covid_data[[score_col_name]][query_cell_matches]
    
    pca_data$score <- scores
    pca_data$original_cell_id <- original_ids
    
    return(pca_data)
}

# _________________
# SingleR - Row 1
# _________________

pca_singler <- projectPCA(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "singler_annotations_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# Create SingleR scores first
singler_scores_matrix <- covid_data$singler_scores
singler_assigned <- covid_data$singler_annotations_merged
singler_conf <- sapply(seq_len(ncol(covid_data)), function(i) {
    cell_type <- singler_assigned[i]
    mean(singler_scores_matrix[i, c("CD14_mono", "CD83_CD14_mono")])
})
names(singler_conf) <- colnames(covid_data)

# Add to PCA data
pca_singler <- add_scores_to_pca(pca_singler, covid_data, normal_data, "singler_scores", "singler_annotations_merged")
pca_singler$dataset <- factor(pca_singler$dataset, levels = c("Query", "Reference"))

# Properly assign scores
query_mask <- pca_singler$dataset == "Query"
query_ids <- pca_singler$original_cell_id[query_mask]
pca_singler$score <- NA
pca_singler$score[query_mask] <- singler_conf[query_ids]

# Add IFN score
ref_cells <- pca_singler$original_cell_id[pca_singler$dataset == "Reference"]
query_cells <- pca_singler$original_cell_id[pca_singler$dataset == "Query"]
full_expr <- cbind(
    assay(normal_data, ASSAY_NAME)[signature_genes_avail, ref_cells, drop = FALSE],
    assay(covid_data, ASSAY_NAME)[signature_genes_avail, query_cells, drop = FALSE]
)
pca_singler$ifn_score <- colMeans(full_expr, na.rm = TRUE)[pca_singler$original_cell_id]

fig_s3_1a <- create_pca_plot_grid(pca_singler, "score", "Delta Score", "SingleR")
fig_s3_1b <- create_pca_plot_grid(pca_singler, "ifn_score", "IFN Signature", "SingleR")
fig_s3_1c <- plot(anomaly_singler,
                  cell_type = CELL_TYPES_TO_PLOT,
                  pc_subset = PC_SUBSET,
                  data_type = "query",
                  n_tree = 500,
                  diagonal_facet = "ridge",
                  upper_facet = "blank",
                  max_cells_query = MAX_CELLS) + 
    ggtitle("SingleR: Anomaly Detection (CD14+ Monocytes)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937"))

# _________________
# Azimuth - Row 2
# _________________

pca_azimuth <- projectPCA(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "azimuth_celltype_l1_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

pca_azimuth$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_azimuth))
pca_azimuth$dataset <- factor(pca_azimuth$dataset, levels = c("Query", "Reference"))

# Add scores
query_mask <- pca_azimuth$dataset == "Query"
query_ids <- pca_azimuth$original_cell_id[query_mask]
pca_azimuth$score <- NA
names(covid_data$azimuth_mapping_score) <- colnames(covid_data)
pca_azimuth$score[query_mask] <- covid_data$azimuth_mapping_score[query_ids]

# Add IFN score
ref_cells <- pca_azimuth$original_cell_id[pca_azimuth$dataset == "Reference"]
query_cells <- pca_azimuth$original_cell_id[pca_azimuth$dataset == "Query"]
full_expr <- cbind(
    assay(normal_data, ASSAY_NAME)[signature_genes_avail, ref_cells, drop = FALSE],
    assay(covid_data, ASSAY_NAME)[signature_genes_avail, query_cells, drop = FALSE]
)
pca_azimuth$ifn_score <- colMeans(full_expr, na.rm = TRUE)[pca_azimuth$original_cell_id]

fig_s3_2a <- create_pca_plot_grid(pca_azimuth, "score", "Prediction Score", "Azimuth")
fig_s3_2b <- create_pca_plot_grid(pca_azimuth, "ifn_score", "IFN Signature", "Azimuth")
fig_s3_2c <- plot(anomaly_azimuth,
                  cell_type = CELL_TYPES_TO_PLOT,
                  pc_subset = PC_SUBSET,
                  data_type = "query",
                  n_tree = 500,
                  diagonal_facet = "ridge",
                  upper_facet = "blank",
                  max_cells_query = MAX_CELLS) + 
    ggtitle("Azimuth: Anomaly Detection (CD14+ Monocytes)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937"))

# ____________________
# CellTypist - Row 3
# ____________________

pca_celltypist <- projectPCA(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "celltypist_predicted_labels_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

pca_celltypist$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_celltypist))
pca_celltypist$dataset <- factor(pca_celltypist$dataset, levels = c("Query", "Reference"))

# Add scores
query_mask <- pca_celltypist$dataset == "Query"
query_ids <- pca_celltypist$original_cell_id[query_mask]
pca_celltypist$score <- NA
names(covid_data$celltypist_conf_score) <- colnames(covid_data)
pca_celltypist$score[query_mask] <- covid_data$celltypist_conf_score[query_ids]

# Add IFN score
ref_cells <- pca_celltypist$original_cell_id[pca_celltypist$dataset == "Reference"]
query_cells <- pca_celltypist$original_cell_id[pca_celltypist$dataset == "Query"]
full_expr <- cbind(
    assay(normal_data, ASSAY_NAME)[signature_genes_avail, ref_cells, drop = FALSE],
    assay(covid_data, ASSAY_NAME)[signature_genes_avail, query_cells, drop = FALSE]
)
pca_celltypist$ifn_score <- colMeans(full_expr, na.rm = TRUE)[pca_celltypist$original_cell_id]

fig_s3_3a <- create_pca_plot_grid(pca_celltypist, "score", "Confidence Score", "CellTypist")
fig_s3_3b <- create_pca_plot_grid(pca_celltypist, "ifn_score", "IFN Signature", "CellTypist")
fig_s3_3c <- plot(anomaly_celltypist,
                  cell_type = CELL_TYPES_TO_PLOT,
                  pc_subset = PC_SUBSET,
                  data_type = "query",
                  n_tree = 500,
                  diagonal_facet = "ridge",
                  upper_facet = "blank",
                  max_cells_query = MAX_CELLS) + 
    ggtitle("CellTypist: Anomaly Detection (CD14+ Monocytes)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937"))

# _________________
# scArches - Row 4
# _________________

pca_scarches <- projectPCA(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "scvi_prediction_merged",
    ref_cell_type_col = "author_cell_type_merged",
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

pca_scarches$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_scarches))
pca_scarches$dataset <- factor(pca_scarches$dataset, levels = c("Query", "Reference"))

# Add scores
query_mask <- pca_scarches$dataset == "Query"
query_ids <- pca_scarches$original_cell_id[query_mask]
pca_scarches$score <- NA
names(covid_data$scvi_confidence) <- colnames(covid_data)
pca_scarches$score[query_mask] <- covid_data$scvi_confidence[query_ids]

# Add IFN score
ref_cells <- pca_scarches$original_cell_id[pca_scarches$dataset == "Reference"]
query_cells <- pca_scarches$original_cell_id[pca_scarches$dataset == "Query"]
full_expr <- cbind(
    assay(normal_data, ASSAY_NAME)[signature_genes_avail, ref_cells, drop = FALSE],
    assay(covid_data, ASSAY_NAME)[signature_genes_avail, query_cells, drop = FALSE]
)
pca_scarches$ifn_score <- colMeans(full_expr, na.rm = TRUE)[pca_scarches$original_cell_id]

fig_s3_4a <- create_pca_plot_grid(pca_scarches, "score", "Uncertainty Score", "scArches")
fig_s3_4b <- create_pca_plot_grid(pca_scarches, "ifn_score", "IFN Signature", "scArches")
fig_s3_4c <- plot(anomaly_scarches,
                  cell_type = CELL_TYPES_TO_PLOT,
                  pc_subset = PC_SUBSET,
                  data_type = "query",
                  n_tree = 500,
                  diagonal_facet = "ridge",
                  upper_facet = "blank",
                  max_cells_query = MAX_CELLS) + 
    ggtitle("scArches: Anomaly Detection (CD14+ Monocytes)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937"))

# _________________
# Combine 4x3 Grid
# _________________

fig_s3_combined <- plot_grid(
    ggmatrix_gtable(fig_s3_1a), ggmatrix_gtable(fig_s3_1b), ggmatrix_gtable(fig_s3_1c),
    ggmatrix_gtable(fig_s3_2a), ggmatrix_gtable(fig_s3_2b), ggmatrix_gtable(fig_s3_2c),
    ggmatrix_gtable(fig_s3_3a), ggmatrix_gtable(fig_s3_3b), ggmatrix_gtable(fig_s3_3c),
    ggmatrix_gtable(fig_s3_4a), ggmatrix_gtable(fig_s3_4b), ggmatrix_gtable(fig_s3_4c),
    nrow = 4, ncol = 3, labels = "AUTO",
    label_size = 12, label_fontface = "bold"
)

ggsave("figures/supp/Fig_S3_pca_comparison.png", fig_s3_combined, width = 18, height = 20, dpi = 600)

# ________
# Summary
# ________

print("Supplementary Figure S3 complete!")
print("Saved in: figures/supp/")

# ____________________________________________
# Table S3: Wilcoxon Test - IFN in Anomalies
# ____________________________________________

# Pre-calculate IFN scores for ALL CD14+ cells in each annotation method

# SingleR
cd14_singler_all <- covid_data[, covid_data$singler_annotations_merged == "CD14 mono"]
ifn_expr_singler <- assay(cd14_singler_all, "logcounts")[signature_genes_avail, ]
ifn_score_singler_all <- colMeans(ifn_expr_singler, na.rm = TRUE)
names(ifn_score_singler_all) <- colnames(cd14_singler_all)

# Azimuth
cd14_azimuth_all <- covid_data[, covid_data$azimuth_celltype_l1_merged == "CD14 mono"]
ifn_expr_azimuth <- assay(cd14_azimuth_all, "logcounts")[signature_genes_avail, ]
ifn_score_azimuth_all <- colMeans(ifn_expr_azimuth, na.rm = TRUE)
names(ifn_score_azimuth_all) <- colnames(cd14_azimuth_all)

# CellTypist
cd14_celltypist_all <- covid_data[, covid_data$celltypist_predicted_labels_merged == "CD14 mono"]
ifn_expr_celltypist <- assay(cd14_celltypist_all, "logcounts")[signature_genes_avail, ]
ifn_score_celltypist_all <- colMeans(ifn_expr_celltypist, na.rm = TRUE)
names(ifn_score_celltypist_all) <- colnames(cd14_celltypist_all)

# scArches
cd14_scarches_all <- covid_data[, covid_data$scvi_prediction_merged == "CD14 mono"]
ifn_expr_scarches <- assay(cd14_scarches_all, "logcounts")[signature_genes_avail, ]
ifn_score_scarches_all <- colMeans(ifn_expr_scarches, na.rm = TRUE)
names(ifn_score_scarches_all) <- colnames(cd14_scarches_all)

# Helper function
calc_ifn_wilcox <- function(method_name, annotation_col, annotation_cell_type, 
                             anomaly_logical, ifn_scores) {
    
    # Get ALL CD14+ cells for this method
    cd14_mask <- covid_data[[annotation_col]] == annotation_cell_type
    cd14_cell_names <- colnames(covid_data)[cd14_mask]
    
    # Get IFN scores for these cells
    ifn_scores_subset <- ifn_scores[cd14_cell_names]
    
    # Get anomaly status
    anomaly_status <- anomaly_logical[cd14_cell_names]
    
    # Split by anomaly
    ifn_anom <- ifn_scores_subset[anomaly_status == TRUE]
    ifn_non_anom <- ifn_scores_subset[anomaly_status == FALSE]
    
    # Remove NAs
    ifn_anom <- ifn_anom[!is.na(ifn_anom)]
    ifn_non_anom <- ifn_non_anom[!is.na(ifn_non_anom)]
    
    if (length(ifn_anom) < 3 || length(ifn_non_anom) < 3) {
        return(NULL)
    }
    
    # Wilcoxon test
    wilcox_result <- wilcox.test(ifn_anom, ifn_non_anom)
    
    p_val_str <- ifelse(wilcox_result$p.value < 0.001, "<0.001",
                       paste0("p=", format(round(wilcox_result$p.value, 4), nsmall = 4)))
    
    row <- data.frame(
        Method = method_name,
        N_Anomalous = length(ifn_anom),
        N_Non_Anomalous = length(ifn_non_anom),
        Median_IFN_Anom = round(median(ifn_anom, na.rm = TRUE), 2),
        Median_IFN_Non_Anom = round(median(ifn_non_anom, na.rm = TRUE), 2),
        Wilcoxon_p_value = p_val_str,
        stringsAsFactors = FALSE
    )
    
    return(row)
}

table_s3_list <- list()

table_s3_list[[1]] <- calc_ifn_wilcox("SingleR", "singler_annotations_merged", "CD14 mono", 
                                       singler_anomaly_logical, ifn_score_singler_all)
table_s3_list[[2]] <- calc_ifn_wilcox("Azimuth", "azimuth_celltype_l1_merged", "CD14 mono", 
                                       azimuth_anomaly_logical, ifn_score_azimuth_all)
table_s3_list[[3]] <- calc_ifn_wilcox("CellTypist", "celltypist_predicted_labels_merged", "CD14 mono", 
                                       celltypist_anomaly_logical, ifn_score_celltypist_all)
table_s3_list[[4]] <- calc_ifn_wilcox("scArches", "scvi_prediction_merged", "CD14 mono", 
                                       scarches_anomaly_logical, ifn_score_scarches_all)

table_s3 <- do.call(rbind, Filter(Negate(is.null), table_s3_list))
rownames(table_s3) <- NULL
print(table_s3)