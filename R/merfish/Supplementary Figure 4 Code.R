# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 5
# -------------------------------------------------------

library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(GGally)
library(ggridges)
library(viridis)
library(cowplot)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\n--- Loading MERFISH Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

set.seed(0)

# _________________
# Fibroblast Setup
# _________________

cat("\n--- Setting up Fibroblast Analysis ---\n")

# For healthy_data: use tier2_merged which has "Fibroblast"
# For dss9_data: use tier2_merged which has "Fibroblast" and "Inflamed Fibroblast"

FIBROBLAST_TYPES <- c("Fibroblast", "Inflamed Fibroblast")

# Count fibroblasts
healthy_fb_count <- sum(healthy_data$tier2_merged == "Fibroblast")
dss9_fb_count <- sum(dss9_data$tier2_merged %in% FIBROBLAST_TYPES)

cat(sprintf("Healthy Fibroblasts: %d\n", healthy_fb_count))
cat(sprintf("DSS9 Fibroblasts: %d\n", dss9_fb_count))

# ___________________________
# Parameters & ECM Signature
# ___________________________

ecm_homeostasis_signature <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")

PC_SUBSET <- 1:3
ASSAY_NAME <- "logcounts"
MAX_CELLS <- 1000
GLOBAL_ECM_LIMIT <- 3.0

# Identify available genes
common_genes <- intersect(rownames(dss9_data), rownames(healthy_data))
signature_genes_avail <- intersect(ecm_homeostasis_signature, common_genes)

cat(sprintf("ECM signature genes available: %s\n", paste(signature_genes_avail, collapse = ", ")))

# ___________________
# Plotting Functions
# ___________________

scatter_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
        geom_point(data = ~subset(., dataset == "Reference"), 
                   aes(color = ecm_score), shape = 1, alpha = 0.5, size = 0.75) +
        geom_point(data = ~subset(., dataset == "Query"), 
                   aes(color = ecm_score), shape = 16, alpha = 0.8, size = 0.75) +
        viridis::scale_color_viridis(
            option = "plasma", 
            direction = 1, 
            limits = c(min(data$ecm_score, na.rm=TRUE), GLOBAL_ECM_LIMIT), 
            oob = scales::squish
        ) + 
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

# _______________________________
# Function to Process One Method
# _______________________________

process_annotation_method <- function(method_name, 
                                      query_annot_col) {
    
    cat(sprintf("\n=== Processing %s ===\n", method_name))
    
    # ===== LEFT: PCA + ECM =====
    cat("Projecting PCA...\n")
    
    pca_output <- scDiagnostics::projectPCA(
        query_data = dss9_data,
        reference_data = healthy_data,
        query_cell_type_col = query_annot_col,
        ref_cell_type_col = "tier2_merged",
        cell_types = "Fibroblast",
        pc_subset = PC_SUBSET,
        assay_name = ASSAY_NAME,
        max_cells_query = MAX_CELLS,
        max_cells_ref = MAX_CELLS
    )
    
    pca_output$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_output))
    pca_output$dataset <- factor(pca_output$dataset, levels = c("Query", "Reference"))
    
    # Calculate ECM Score
    cat("Calculating ECM score...\n")
    
    ref_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Reference"]
    query_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Query"]
    
    full_expr_matrix <- cbind(
        assay(healthy_data, ASSAY_NAME)[signature_genes_avail, ref_cells_toplot, drop = FALSE],
        assay(dss9_data, ASSAY_NAME)[signature_genes_avail, query_cells_toplot, drop = FALSE]
    )
    
    if(length(signature_genes_avail) > 0) {
        raw_score <- colMeans(full_expr_matrix, na.rm = TRUE)
        capped_score <- pmin(raw_score, 1.5)
        final_score <- capped_score * 2
        pca_output$ecm_score <- final_score[match(pca_output$original_cell_id, names(final_score))]
    } else {
        pca_output$ecm_score <- 0
    }
    
    # Create legend
    legend_plot <- ggplot(pca_output, aes(x = PC1, y = PC2, color = ecm_score)) +
        geom_point() +
        viridis::scale_color_viridis(
            option = "plasma", 
            direction = 1, 
            name = "ECM\nScore",
            limits = c(min(pca_output$ecm_score, na.rm=TRUE), GLOBAL_ECM_LIMIT),
            oob = scales::squish
        ) +
        theme(legend.position = "right", legend.box = "vertical", 
              legend.key = element_rect(fill = NA, color = NA))
    
    plot_legend <- GGally::grab_legend(legend_plot)
    
    # Create PCA plot
    pc_plot_names <- paste0(
        "PC", PC_SUBSET, " (",
        sprintf("%.1f%%", attributes(reducedDim(healthy_data, "PCA"))[["percentVar"]][PC_SUBSET]), ")"
    )
    
    cat("Creating PCA+ECM plot...\n")
    
    pca_ecm_plot <- GGally::ggpairs(
        pca_output,
        columns = paste0("PC", PC_SUBSET),
        columnLabels = pc_plot_names,
        mapping = aes(color = ecm_score, shape = dataset), 
        lower = list(continuous = scatter_fn),
        diag = list(continuous = ridge_fn),
        upper = list(continuous = blank_fn),
        progress = FALSE,
        legend = plot_legend,
        title = paste0(method_name, " - PCA + ECM")
    ) + 
    theme(
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black", face = "bold", size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937")
    )
    
    # ===== RIGHT: Anomaly Detection =====
    cat("Running anomaly detection...\n")
    
    anomaly_output <- detectAnomaly(
        query_data = dss9_data,
        reference_data = healthy_data, 
        query_cell_type_col = query_annot_col, 
        ref_cell_type_col = "tier2_merged", 
        cell_types = "Fibroblast", 
        pc_subset = PC_SUBSET,
        n_tree = 500,
        anomaly_threshold = 0.5, 
        assay_name = ASSAY_NAME,
        max_cells_query = NULL, 
        max_cells_ref = NULL
    )
    
    cat("Creating anomaly detection plot...\n")
    
    anomaly_plot <- plot(
        anomaly_output,
        cell_type = "Fibroblast", 
        pc_subset = PC_SUBSET,
        data_type = "query",  
        n_tree = 500,
        diagonal_facet = "ridge",
        upper_facet = "blank",
        max_cells_query = 2000 
    ) + 
        ggtitle(paste0(method_name, " - Anomaly Detection")) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1F2937"),
            strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            strip.text = element_text(color = "black", face = "bold", size = 9)
        )
    
    return(list(pca_plot = pca_ecm_plot, anomaly_plot = anomaly_plot))
}

# ____________________
# Process All Methods
# ____________________

result_singler <- process_annotation_method("SingleR", "singler_annotations_merged")
result_azimuth <- process_annotation_method("Azimuth", "azimuth_celltype_l1_merged")
result_celltypist <- process_annotation_method("CellTypist", "celltypist_predicted_labels_merged")
result_scarches <- process_annotation_method("scArches", "scvi_prediction_merged")

# ______________________
# Combine into 4x2 Grid
# ______________________

fig_s4_combined <- plot_grid(
    ggmatrix_gtable(result_singler$pca_plot), ggmatrix_gtable(result_singler$anomaly_plot),
    ggmatrix_gtable(result_azimuth$pca_plot), ggmatrix_gtable(result_azimuth$anomaly_plot),
    ggmatrix_gtable(result_celltypist$pca_plot), ggmatrix_gtable(result_celltypist$anomaly_plot),
    ggmatrix_gtable(result_scarches$pca_plot), ggmatrix_gtable(result_scarches$anomaly_plot),
    nrow = 4, ncol = 2, labels = "AUTO",
    label_size = 12, label_fontface = "bold",
    rel_widths = c(1, 1)
)

# _________________
# Display and Save
# _________________

cat("Creating output directory...\n")
dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)

cat("Saving combined figure...\n")
ggsave("figures/supp/merfish/Fig_S4_pca_ecm_anomaly.png", 
       fig_s4_combined, 
       width = 18, 
       height = 20, 
       dpi = 600)

cat("✓ Figure saved\n")
