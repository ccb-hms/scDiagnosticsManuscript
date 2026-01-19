# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 5 Code
# -------------------------------------------------------

library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(scater)
library(BiocSingular)
library(cowplot)
library(dplyr)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\n--- Loading MERFISH Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

set.seed(0)

# _____________________
# ECM Signature Setup
# _____________________

cat("\n--- Setting up ECM Gene Analysis ---\n")

ecm_genes <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")

cat(sprintf("ECM genes to analyze: %s\n", paste(ecm_genes, collapse = ", ")))

# Check availability
common_genes <- intersect(rownames(dss9_data), rownames(healthy_data))
available_ecm_genes <- intersect(ecm_genes, common_genes)

cat(sprintf("Available ECM genes: %s\n", paste(available_ecm_genes, collapse = ", ")))

cat("\n--- Subsetting to ECM genes only ---\n")

healthy_data_ecm <- healthy_data[available_ecm_genes, ]
dss9_data_ecm <- dss9_data[available_ecm_genes, ]

cat(sprintf("Healthy data dims: %d genes x %d cells\n", nrow(healthy_data_ecm), ncol(healthy_data_ecm)))
cat(sprintf("DSS9 data dims: %d genes x %d cells\n", nrow(dss9_data_ecm), ncol(dss9_data_ecm)))

cat("\n--- Computing PCA on reference (all ECM genes, no HVG filtering) ---\n")

healthy_data_ecm <- scater::runPCA(
    healthy_data_ecm,
    subset_row = rownames(healthy_data_ecm),  # Use ALL genes, not HVGs
    ncomponents = 3,
    BSPARAM = BiocSingular::IrlbaParam(),
    name = "PCA"
)

cat("✓ PCA computed on reference with all ECM genes\n")
cat(sprintf("PCA dimensions: %s\n", paste(dim(reducedDim(healthy_data_ecm, "PCA")), collapse = " x ")))

# _________________________________________
# Function to Reorder Genes Alphabetically
# _________________________________________

reorder_genes_alphabetically <- function(plot_obj) {
    # Extract plot data
    plot_data <- plot_obj$data
    
    # Get unique genes in alphabetical order (reversed for y-axis)
    genes_alpha <- sort(unique(as.character(plot_data$gene)))
    genes_alpha_reversed <- rev(genes_alpha)
    
    # Update the factor levels
    plot_obj$data$gene <- factor(plot_obj$data$gene, levels = genes_alpha_reversed)
    
    cat(sprintf("Gene order (alphabetical): %s\n", paste(genes_alpha, collapse = ", ")))
    
    return(plot_obj)
}

# ________________________________
# Function to Process One Method
# ________________________________

process_annotation_method <- function(method_name, 
                                      query_annot_col) {
    
    cat(sprintf("\n=== Processing %s ===\n", method_name))
    
    # ===== Calculate Gene Shifts =====
    cat("Calculating gene shifts with anomaly detection...\n")
    
    gene_shifts <- calculateGeneShifts(
        query_data = dss9_data_ecm,
        reference_data = healthy_data_ecm,
        query_cell_type_col = query_annot_col,
        ref_cell_type_col = "tier2_merged",
        cell_types = "Fibroblast",
        pc_subset = 1:3,
        n_top_loadings = 5,
        p_value_threshold = 0.05,
        adjust_method = "fdr",
        assay_name = "logcounts",
        detect_anomalies = TRUE,
        anomaly_comparison = TRUE,
        anomaly_threshold = 0.5,
        n_tree = 500,
        max_cells_query = NULL,
        max_cells_ref = NULL
    )
    
    # Verify anomaly detection worked
    has_anomaly <- "anomaly_results" %in% names(gene_shifts)
    cat(sprintf("Anomaly detection: %s\n", ifelse(has_anomaly, "✓ Success", "✗ Failed")))
    
    # ===== CREATE BARPLOT =====
    cat("Creating barplot...\n")
    
    barplot_plot <- plot(
        x = gene_shifts,
        cell_type = "Fibroblast",
        pc_subset = 1:3,
        plot_type = "barplot",
        plot_by = "p_adjusted",
        n_genes = length(available_ecm_genes),
        significance_threshold = 0.05,
        show_anomalies = TRUE,
        pseudo_bulk = TRUE,
        show_all_query = TRUE
    ) +
        ggtitle(paste0(method_name, " - Fibroblasts (ECM Genes)")) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1F2937"),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.4, "cm")
        )
    
    # ===== REORDER GENES ALPHABETICALLY =====
    cat("Reordering genes alphabetically...\n")
    barplot_plot <- reorder_genes_alphabetically(barplot_plot)
    
    return(barplot_plot)
}

# ____________________
# Process All Methods
# ____________________

result_azimuth <- process_annotation_method("Azimuth", "azimuth_celltype_l1_merged")
result_singler <- process_annotation_method("SingleR", "singler_annotations_merged")
result_celltypist <- process_annotation_method("CellTypist", "celltypist_predicted_labels_merged")
result_scarches <- process_annotation_method("scArches", "scvi_prediction_merged")

# ______________________
# Combine into 2x2 Grid
# ______________________

fig_s5_combined <- plot_grid(
    result_azimuth, result_singler,
    result_celltypist, result_scarches,
    nrow = 2, ncol = 2, labels = "AUTO",
    label_size = 12, label_fontface = "bold"
)

# _________________
# Display and Save
# _________________

cat("Creating output directory...\n")
dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)

cat("Saving combined figure...\n")
ggsave("figures/supp/merfish/Fig_S5_ecm_barplot.png", 
       fig_s5_combined, 
       width = 14, 
       height = 10, 
       dpi = 600)
cat("✓ Figure saved\n")

# _____________________________________________________________
# Supplementary Table S5: Col1a2 Expression by Anomaly Status
# _____________________________________________________________

cat("\n--- Generating Supplementary Table 5: Col1a2 Anomaly Analysis ---\n")

# Function to detect anomalies and compare Col1a2 expression
analyze_col1a2_by_anomaly <- function(query_data, query_annot_col, 
                                      reference_data, ref_annot_col,
                                      method_name) {
    
    cat(sprintf("\nProcessing %s...\n", method_name))
    
    # Detect anomalies using detectAnomaly function
    anomaly_output <- detectAnomaly(
        reference_data = reference_data,
        query_data = query_data,
        ref_cell_type_col = ref_annot_col,
        query_cell_type_col = query_annot_col,
        cell_types = "Fibroblast",
        pc_subset = 1:3,
        n_tree = 500,
        anomaly_threshold = 0.5,
        assay_name = "logcounts", 
        max_cells_ref = NULL, 
        max_cells_query = NULL
    )
    
    # Extract anomaly status for Fibroblast cell type
    anomaly_status <- anomaly_output[["Fibroblast"]][["query_anomaly"]]
    
    # Get Col1a2 expression for all cells
    col1a2_expr <- assay(query_data, "logcounts")["Col1a2", ]
    
    # Get Fibroblast indices in query data
    fib_idx <- colData(query_data)[[query_annot_col]] == "Fibroblast"
    col1a2_fib <- col1a2_expr[fib_idx]
    
    # Separate by anomaly status
    anom_cells <- col1a2_fib[anomaly_status]
    non_anom_cells <- col1a2_fib[!anomaly_status]
    
    # Calculate stats
    n_anom <- length(anom_cells)
    n_non_anom <- length(non_anom_cells)
    mean_anom <- mean(anom_cells)
    mean_non_anom <- mean(non_anom_cells)
    sd_anom <- sd(anom_cells)
    sd_non_anom <- sd(non_anom_cells)
    
    # Perform t-test
    t_test <- t.test(anom_cells, non_anom_cells)
    
    return(data.frame(
        Method = method_name,
        N_Anomalous = n_anom,
        N_Normal = n_non_anom,
        Mean_Anomalous = round(mean_anom, 3),
        SD_Anomalous = round(sd_anom, 3),
        Mean_Normal = round(mean_non_anom, 3),
        SD_Normal = round(sd_non_anom, 3),
        Mean_Diff = round(mean_anom - mean_non_anom, 3),
        T_Statistic = round(t_test$statistic, 3),
        P_Value = format(t_test$p.value, scientific = TRUE, digits = 3),
        Significant = ifelse(t_test$p.value < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
    ))
}

# Calculate for all methods
table_azimuth <- analyze_col1a2_by_anomaly(dss9_data_ecm, "azimuth_celltype_l1_merged", 
                                            healthy_data_ecm, "tier2_merged", "Azimuth")
table_singler <- analyze_col1a2_by_anomaly(dss9_data_ecm, "singler_annotations_merged", 
                                            healthy_data_ecm, "tier2_merged", "SingleR")
table_celltypist <- analyze_col1a2_by_anomaly(dss9_data_ecm, "celltypist_predicted_labels_merged", 
                                               healthy_data_ecm, "tier2_merged", "CellTypist")
table_scarches <- analyze_col1a2_by_anomaly(dss9_data_ecm, "scvi_prediction_merged", 
                                             healthy_data_ecm, "tier2_merged", "scArches")

# Combine into table
supp_table_5 <- rbind(table_azimuth, table_singler, table_celltypist, table_scarches)

cat("\n\nSupplementary Table 5: Col1a2 Expression - Anomalous vs. Normal Cells\n")
print(supp_table_5)