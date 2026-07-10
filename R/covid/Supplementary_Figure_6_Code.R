# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 6 Code
# -----------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scDiagnostics)

# -----------------------------------------------

# __________
# Load data
# __________

covid_data <- readRDS("data/covid/covid_data_sce.rds")
normal_data <- readRDS("data/covid/normal_data_sce.rds")

set.seed(0)

# ________________________
# Calculate IFN Signature
# ________________________

yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

# ____________
# Setup theme
# ____________

theme_set(theme_minimal(base_size = 20) + theme(
  axis.title = element_text(size = 20, face = "bold"),
  axis.text = element_text(size = 16),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 20, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 18, face = "bold")
))

# ________________________________________________________
# Function: Extract gene order from scDiagnostics results
# ________________________________________________________

extract_geneshifts_order <- function(gene_shifts_obj, pc_subset = 1:5) {
    
    # Combine results from all PCs and sort by p-value
    all_results <- data.frame()
    
    for (pc in pc_subset) {
        pc_name <- paste0("PC", pc)
        if (pc_name %in% names(gene_shifts_obj)) {
            pc_results <- gene_shifts_obj[[pc_name]]
            # Check if it's a data frame with rows
            if (is.data.frame(pc_results) && length(pc_results) > 0) {
                if (!is.null(pc_results) && nrow(pc_results) > 0) {
                    all_results <- rbind(all_results, pc_results)
                }
            }
        }
    }
    
    if (nrow(all_results) == 0) {
        warning("No results found in gene shifts object")
        return(NULL)
    }
    
    # Get unique genes sorted by p-value (ascending)
    all_results <- all_results[order(all_results$p_adjusted), ]
    unique_genes <- unique(all_results$gene)
    
    # Return top 25 genes
    return(head(unique_genes, 25))
}


# __________________________________________
# Function: Create barplot from gene shifts 
# __________________________________________

create_geneshifts_barplot <- function(gene_shifts_obj, cell_type, genes_order, method_name) {
    
    # Get pseudo-bulk expression data
    cell_meta <- gene_shifts_obj$cell_metadata
    expr_data <- gene_shifts_obj$expression_data
    
    # Separate by dataset and anomaly status
    ref_cells <- cell_meta$cell_id[cell_meta$dataset == "Reference" & 
                                    cell_meta$anomaly_status == "Normal"]
    query_normal_cells <- cell_meta$cell_id[cell_meta$dataset == "Query" & 
                                             cell_meta$anomaly_status == "Normal"]
    query_anom_cells <- cell_meta$cell_id[cell_meta$dataset == "Query" & 
                                          cell_meta$anomaly_status == "Anomaly"]
    query_all_cells <- cell_meta$cell_id[cell_meta$dataset == "Query"]
    
    # Calculate pseudo-bulk means
    ref_mean <- rowMeans(expr_data[, ref_cells], na.rm = TRUE)
    query_normal_mean <- rowMeans(expr_data[, query_normal_cells], na.rm = TRUE)
    query_anom_mean <- rowMeans(expr_data[, query_anom_cells], na.rm = TRUE)
    query_all_mean <- rowMeans(expr_data[, query_all_cells], na.rm = TRUE)
    
    # Calculate fold changes
    lfc_normal <- query_normal_mean - ref_mean
    lfc_anom <- query_anom_mean - ref_mean
    lfc_all <- query_all_mean - ref_mean
    
    # Filter to genes in order
    genes_avail <- genes_order[genes_order %in% names(lfc_normal)]
    
    # Create data frame
    plot_data <- data.frame(
      gene = rep(genes_avail, 3),
      log2fc = c(lfc_all[genes_avail], lfc_normal[genes_avail], lfc_anom[genes_avail]),
      group = factor(rep(c("all_query", "normal", "anomaly"), each = length(genes_avail)),
                     levels = c("anomaly", "normal", "all_query")),
      stringsAsFactors = FALSE
    )
    
    # Plot with FIXED gene order
    p <- ggplot(plot_data, aes(x = log2fc, y = factor(gene, levels = rev(genes_avail)), fill = group)) +
        geom_col(position = position_dodge(width = 0.7), alpha = 0.8, width = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
        scale_fill_manual(
            values = c(
                "all_query" = "#9E9E9E",
                "normal" = "#2ecc71",
                "anomaly" = "#e74c3c"
            ),
            labels = c(
                "all_query" = "All Query vs Ref",
                "normal" = "Query Non-Anomaly vs Ref",
                "anomaly" = "Query Anomaly vs Ref"
            ),
            breaks = c("all_query", "normal", "anomaly"),
            name = "Comparison"
        ) +
        xlab("Log2 Fold Change (Pseudo-Bulk)") +
        ylab("") +
        ggtitle(paste(method_name, "- scDiagnostics Anomaly Detection")) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "#1F2937", family = "sans"),
            axis.title = element_text(size = 20, face = "bold", color = "#374151", family = "sans"),
            axis.text.y = element_text(size = 16, color = "#4B5563", family = "sans", face = "bold"),
            axis.text.x = element_text(size = 16, color = "#4B5563", family = "sans"),
            legend.position = "none", # Removed to extract to bottom later
            panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.2),
            panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
            plot.margin = margin(t = 35, r = 10, b = 10, l = 0, "pt"), # Reduced left margin, big top margin
            panel.background = element_rect(fill = "#FAFBFC", color = NA)
        )
    
    return(p)
}

# ____________________________________
# Function: Create gene quantile plot 
# ____________________________________

create_quantile_gene_plot <- function(covid_cd14, normal_cd14, genes, method_name, confidence_col, 
                                      cell_type_for_scores = NULL, score_direction = "higher_is_better",
                                      gene_order = NULL,
                                      assay_name = "logcounts") {
    
    # Get confidence scores for query cells
    if (!is.null(cell_type_for_scores)) {
        query_scores <- as.numeric(covid_cd14[[confidence_col]][, cell_type_for_scores])
    } else {
        query_scores <- as.numeric(covid_cd14[[confidence_col]])
    }
    
    # Remove NA values
    valid_idx <- !is.na(query_scores)
    query_scores_valid <- query_scores[valid_idx]
    query_cell_names_valid <- colnames(covid_cd14)[valid_idx]
    
    # Define high/low confidence based on score direction
    if (score_direction == "higher_is_better") {
        q25 <- quantile(query_scores_valid, 0.25, na.rm = TRUE)
        q75 <- quantile(query_scores_valid, 0.75, na.rm = TRUE)
        low_conf_idx <- query_scores_valid <= q25
        high_conf_idx <- query_scores_valid >= q75
    } else if (score_direction == "lower_is_better") {
        q25 <- quantile(query_scores_valid, 0.25, na.rm = TRUE)
        q75 <- quantile(query_scores_valid, 0.75, na.rm = TRUE)
        low_conf_idx <- query_scores_valid >= q75
        high_conf_idx <- query_scores_valid <= q25
    }
    
    low_conf_cells <- query_cell_names_valid[low_conf_idx]
    high_conf_cells <- query_cell_names_valid[high_conf_idx]
    all_query_cells <- query_cell_names_valid
    
    # Filter genes that exist
    genes_available <- genes[genes %in% rownames(covid_cd14) & genes %in% rownames(normal_cd14)]
    
    # Calculate mean expression for each group and each gene
    ref_mean <- rowMeans(assay(normal_cd14[genes_available, ], assay_name), na.rm = TRUE)
    all_query_mean <- rowMeans(assay(covid_cd14[genes_available, all_query_cells], assay_name), na.rm = TRUE)
    low_conf_mean <- rowMeans(assay(covid_cd14[genes_available, low_conf_cells], assay_name), na.rm = TRUE)
    high_conf_mean <- rowMeans(assay(covid_cd14[genes_available, high_conf_cells], assay_name), na.rm = TRUE)
    
    # Calculate log2 fold changes
    lfc_all <- all_query_mean - ref_mean
    lfc_low <- low_conf_mean - ref_mean
    lfc_high <- high_conf_mean - ref_mean
    
    # Use simple keys without line breaks
    plot_data <- data.frame(
      gene = rep(genes_available, 3),
      log2fc = c(lfc_all, lfc_low, lfc_high),
      group = factor(rep(c("all_query", "low_conf", "high_conf"), 
                        each = length(genes_available)),
                    levels = c("high_conf", "low_conf", "all_query")),
      stringsAsFactors = FALSE
    )
    
    # Create display labels with line breaks
    if (score_direction == "higher_is_better") {
        low_label <- "Query Low Confidence\n(≤25th %ile) vs Ref"
        high_label <- "Query High Confidence\n(≥75th %ile) vs Ref"
    } else {
        low_label <- "Low Confidence\n(≥75th %ile) vs Ref"
        high_label <- "High Confidence\n(≤25th %ile) vs Ref"
    }
    
    # Use gene_order if provided, otherwise order by lfc
    if (!is.null(gene_order)) {
        genes_in_order <- intersect(gene_order, genes_available)
        y_order <- rev(genes_in_order)
        plot_data$gene <- factor(plot_data$gene, levels = genes_in_order)
    } else {
        y_order <- NULL  
    }
    
    # Plot
    if (!is.null(gene_order)) {
        p <- ggplot(plot_data, aes(x = log2fc, y = gene, fill = group))
    } else {
        p <- ggplot(plot_data, aes(x = log2fc, y = reorder(gene, log2fc), fill = group))
    }
    
    p <- p +
        geom_col(position = position_dodge(width = 0.7), alpha = 0.8, width = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
        scale_fill_manual(
            values = c(
                "all_query" = "#5B7C99",
                "low_conf" = "#C1666B",
                "high_conf" = "#52A552"
            ),
            labels = c(
                "all_query" = "All Query vs Ref",
                "low_conf" = low_label,
                "high_conf" = high_label
            ),
            breaks = c("all_query", "low_conf", "high_conf"),
            name = "Comparison"
        ) +
        xlab("Log2 Fold Change (Pseudo-Bulk)") +
        ylab("") +
        ggtitle(paste(method_name, "- Annotation Tool Scores")) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "#1F2937", family = "sans"),
            axis.title = element_text(size = 20, face = "bold", color = "#374151", family = "sans"),
            axis.text.y = element_text(size = 16, color = "#4B5563", family = "sans", face = "bold"),
            axis.text.x = element_text(size = 16, color = "#4B5563", family = "sans"),
            legend.position = "none", # Removed to extract to bottom later
            panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.2),
            panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
            plot.margin = margin(t = 35, r = 5, b = 10, l = 10, "pt"), # Reduced right margin, big top margin
            panel.background = element_rect(fill = "#FAFBFC", color = NA)
        )
    
    return(p)
}

# _________________
# Azimuth - Row 1
# _________________

cd14_azimuth <- covid_data[, covid_data$azimuth_celltype_l1 == "CD14_mono"]
normal_cd14 <- normal_data[, normal_data$author_cell_type == "CD14_mono"]

gene_shifts_azimuth <- calculateGeneShifts(
    query_data = covid_data[yoshida_ifn_signature,],
    reference_data = normal_data[yoshida_ifn_signature,],
    query_cell_type_col = "azimuth_celltype_l1",
    ref_cell_type_col = "author_cell_type",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    p_value_threshold = 0,
    genes_to_analyze = yoshida_ifn_signature,
    detect_anomalies = TRUE,
    anomaly_comparison = TRUE,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

gene_order_azimuth <- extract_geneshifts_order(gene_shifts_azimuth, pc_subset = 1:5)

fig_s6_1a <- create_quantile_gene_plot(cd14_azimuth, normal_cd14, yoshida_ifn_signature, 
                                       "Azimuth", "azimuth_mapping_score",
                                       score_direction = "higher_is_better",
                                       gene_order = NULL,
                                       assay_name = "logcounts")

fig_s6_1b <- create_geneshifts_barplot(gene_shifts_azimuth, "CD14_mono", gene_order_azimuth, "Azimuth")

# _________________
# SingleR - Row 2
# _________________

cd14_singler <- covid_data[, covid_data$singler_annotations_merged == "CD14 mono"]

gene_shifts_singler <- calculateGeneShifts(
    query_data = covid_data[yoshida_ifn_signature,],
    reference_data = normal_data[yoshida_ifn_signature,],
    query_cell_type_col = "singler_annotations",
    ref_cell_type_col = "author_cell_type",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    p_value_threshold = 0,
    detect_anomalies = TRUE,
    anomaly_comparison = TRUE,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

gene_order_singler <- extract_geneshifts_order(gene_shifts_singler, pc_subset = 1:5)

fig_s6_2a <- create_quantile_gene_plot(cd14_singler, normal_cd14, yoshida_ifn_signature, 
                                        "SingleR", "singler_scores", 
                                        cell_type_for_scores = "CD14_mono",
                                        score_direction = "higher_is_better",
                                        gene_order = rev(gene_order_singler))

fig_s6_2b <- create_geneshifts_barplot(gene_shifts_singler, "CD14_mono", gene_order_singler, "SingleR")

# ____________________
# CellTypist - Row 3
# ____________________

cd14_celltypist <- covid_data[, covid_data$celltypist_predicted_labels_merged == "CD14 mono"]

gene_shifts_celltypist <- calculateGeneShifts(
    query_data = covid_data[yoshida_ifn_signature,],
    reference_data = normal_data[yoshida_ifn_signature,],
    query_cell_type_col = "celltypist_predicted_labels",
    ref_cell_type_col = "author_cell_type",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    p_value_threshold = 0,
    genes_to_analyze = yoshida_ifn_signature,
    detect_anomalies = TRUE,
    anomaly_comparison = TRUE,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = 5000,
    max_cells_ref = 5000
)

gene_order_celltypist <- extract_geneshifts_order(gene_shifts_celltypist, pc_subset = 1:5)

fig_s6_3a <- create_quantile_gene_plot(cd14_celltypist, normal_cd14, yoshida_ifn_signature, 
                                        "CellTypist", "celltypist_conf_score",
                                        score_direction = "higher_is_better",
                                        gene_order = rev(gene_order_celltypist))

fig_s6_3b <- create_geneshifts_barplot(gene_shifts_celltypist, "CD14_mono", gene_order_celltypist, "CellTypist")

# _________________
# scArches - Row 4
# _________________

cd14_scarches <- covid_data[, covid_data$scvi_prediction_merged == "CD14 mono"]

gene_shifts_scarches <- calculateGeneShifts(
    query_data = covid_data[yoshida_ifn_signature,],
    reference_data = normal_data[yoshida_ifn_signature,],
    query_cell_type_col = "scvi_prediction",
    ref_cell_type_col = "author_cell_type",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    p_value_threshold = 0,
    genes_to_analyze = yoshida_ifn_signature,
    detect_anomalies = TRUE,
    anomaly_comparison = TRUE,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = 5000,
    max_cells_ref = 5000
)

gene_order_scarches <- extract_geneshifts_order(gene_shifts_scarches, pc_subset = 1:5)

fig_s6_4a <- create_quantile_gene_plot(cd14_scarches, normal_cd14, yoshida_ifn_signature, 
                                        "scArches", "scvi_confidence",
                                        score_direction = "lower_is_better",
                                        gene_order = rev(gene_order_scarches))

fig_s6_4b <- create_geneshifts_barplot(gene_shifts_scarches, "CD14_mono", gene_order_scarches, "scArches")

# _________________
# Combine 4x2 grid
# _________________

# Build the main grid (legends are stripped out via the updated functions)
grid_plots <- plot_grid(
    fig_s6_1a, fig_s6_1b,
    fig_s6_2a, fig_s6_2b,
    fig_s6_3a, fig_s6_3b,
    fig_s6_4a, fig_s6_4b,
    nrow = 4, ncol = 2, labels = "AUTO",
    label_size = 26, label_fontface = "bold", 
    label_x = 0, label_y = 1, # Snaps labels safely into the top margin gap
    rel_widths = c(1, 1)
)

# Extract legends for the bottom
legend_left <- cowplot::get_legend(
    fig_s6_1a + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 18)) +
    guides(fill = guide_legend(nrow = 2, reverse = TRUE))
)
legend_right <- cowplot::get_legend(
    fig_s6_1b + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 18)) +
    guides(fill = guide_legend(nrow = 2, reverse = TRUE))
)

# Place legends side-by-side
bottom_legends <- plot_grid(legend_left, legend_right, ncol = 2)

# Combine the grid and the legends vertically
fig_s6_combined <- plot_grid(grid_plots, bottom_legends, ncol = 1, rel_heights = c(1, 0.08))

# Width set to 20 instead of 28 to pull columns closer together
ggsave("figures/supp/covid/Fig_S6_gene_expression_shifts.png", fig_s6_combined, width = 20, height = 30, dpi = 600)

# ________
# Summary
# ________

print("Supplementary Figure S6 complete!")
print("Saved in: figures/supp/covid/")

# _____________________________________________________________________
# Table S6: Correlation - Annotation Quantile vs scDiagnostics Anomaly
# _____________________________________________________________________

create_correlation_row <- function(gene_order, cd14_data, normal_cd14_data, 
                                   gene_shifts_obj, method_name, 
                                   confidence_col, cell_type_for_scores = NULL,
                                   score_direction = "higher_is_better") {
    
    # Get genes
    genes_avail <- gene_order[gene_order %in% rownames(cd14_data)]
    
    # ===== LEFT PLOT: Annotation Quantile =====
    
    # Get confidence scores
    if (!is.null(cell_type_for_scores)) {
        query_scores <- as.numeric(cd14_data[[confidence_col]][, cell_type_for_scores])
    } else {
        query_scores <- as.numeric(cd14_data[[confidence_col]])
    }
    
    valid_idx <- !is.na(query_scores)
    query_scores_valid <- query_scores[valid_idx]
    query_cell_names_valid <- colnames(cd14_data)[valid_idx]
    
    # Define high/low
    if (score_direction == "higher_is_better") {
        q25 <- quantile(query_scores_valid, 0.25, na.rm = TRUE)
        q75 <- quantile(query_scores_valid, 0.75, na.rm = TRUE)
        low_conf_idx <- query_scores_valid <= q25
        high_conf_idx <- query_scores_valid >= q75
    } else {
        q25 <- quantile(query_scores_valid, 0.25, na.rm = TRUE)
        q75 <- quantile(query_scores_valid, 0.75, na.rm = TRUE)
        low_conf_idx <- query_scores_valid >= q75
        high_conf_idx <- query_scores_valid <= q25
    }
    
    low_conf_cells <- query_cell_names_valid[low_conf_idx]
    high_conf_cells <- query_cell_names_valid[high_conf_idx]
    
    # Calculate fold changes for annotation quantile
    ref_mean <- rowMeans(assay(normal_cd14_data[genes_avail, ]), na.rm = TRUE)
    low_conf_mean <- rowMeans(assay(cd14_data[genes_avail, low_conf_cells]), na.rm = TRUE)
    high_conf_mean <- rowMeans(assay(cd14_data[genes_avail, high_conf_cells]), na.rm = TRUE)
    
    lfc_annot <- high_conf_mean - low_conf_mean
    
    # ===== RIGHT PLOT: scDiagnostics Anomaly =====
    
    cell_meta <- gene_shifts_obj$cell_metadata
    expr_data <- gene_shifts_obj$expression_data
    
    ref_cells <- cell_meta$cell_id[cell_meta$dataset == "Reference" & 
                                    cell_meta$anomaly_status == "Normal"]
    query_normal_cells <- cell_meta$cell_id[cell_meta$dataset == "Query" & 
                                             cell_meta$anomaly_status == "Normal"]
    query_anom_cells <- cell_meta$cell_id[cell_meta$dataset == "Query" & 
                                          cell_meta$anomaly_status == "Anomaly"]
    
    ref_mean_anom <- rowMeans(expr_data[, ref_cells], na.rm = TRUE)
    query_normal_mean_anom <- rowMeans(expr_data[, query_normal_cells], na.rm = TRUE)
    query_anom_mean_anom <- rowMeans(expr_data[, query_anom_cells], na.rm = TRUE)
    
    lfc_anom_vs_normal <- query_anom_mean_anom - query_normal_mean_anom
    
    # ===== CORRELATION =====
    
    # Match genes
    genes_match <- intersect(names(lfc_annot), names(lfc_anom_vs_normal))
    
    if (length(genes_match) < 3) {
        warning(paste(method_name, ": Not enough matching genes for correlation"))
        return(NULL)
    }
    
    lfc_annot_matched <- lfc_annot[genes_match]
    lfc_anom_matched <- lfc_anom_vs_normal[genes_match]
    
    # Spearman correlation
    cor_result <- cor.test(lfc_annot_matched, lfc_anom_matched, method = "spearman")
    
    p_val_str <- ifelse(cor_result$p.value < 0.001, "<0.001",
                       paste0("p=", format(round(cor_result$p.value, 4), nsmall = 4)))
    
    row <- data.frame(
        Method = method_name,
        N_Genes = length(genes_match),
        Spearman_r = round(cor_result$estimate, 3),
        p_value = p_val_str,
        stringsAsFactors = FALSE
    )
    
    return(row)
}

# Generate table S6
table_s6_list <- list(
    create_correlation_row(gene_order_azimuth, cd14_azimuth, normal_cd14,
                          gene_shifts_azimuth, "Azimuth", "azimuth_mapping_score",
                          score_direction = "higher_is_better"),
    
    create_correlation_row(gene_order_singler, cd14_singler, normal_cd14,
                          gene_shifts_singler, "SingleR", "singler_scores",
                          cell_type_for_scores = "CD14_mono",
                          score_direction = "higher_is_better"),
    
    create_correlation_row(gene_order_celltypist, cd14_celltypist, normal_cd14,
                          gene_shifts_celltypist, "CellTypist", "celltypist_conf_score",
                          score_direction = "higher_is_better"),
    
    create_correlation_row(gene_order_scarches, cd14_scarches, normal_cd14,
                          gene_shifts_scarches, "scArches", "scvi_confidence",
                          score_direction = "lower_is_better")
)

table_s6 <- do.call(rbind, Filter(Negate(is.null), table_s6_list))
rownames(table_s6) <- NULL
print(table_s6)