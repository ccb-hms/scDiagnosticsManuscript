# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 6 Code
# -------------------------------------------------------

library(SingleCellExperiment)
library(scDiagnostics)
library(scran)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(cowplot)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\n--- Loading MERFISH Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

set.seed(0)

# _______________________________________________________
# Function to Create Heatmaps for Each Annotation Method
# _______________________________________________________

create_annotation_heatmaps <- function(method_name, query_annot_col) {
    
    cat(sprintf("\n=== Processing %s ===\n", method_name))
    
    # _______________________________
    # Get cell groups for comparison
    # _______________________________
    
    cat("\nIdentifying cell groups...\n")
    
    # 1. Inflamed Fibro annotated as SMC
    inflamed_to_smc <- dss9_data[, 
        dss9_data$tier2_merged == "Inflamed Fibroblast" & 
        dss9_data[[query_annot_col]] == "Smooth Muscle"
    ]
    
    # 2. Inflamed Fibro annotated as Fibro
    inflamed_to_fibro <- dss9_data[, 
        dss9_data$tier2_merged == "Inflamed Fibroblast" & 
        dss9_data[[query_annot_col]] == "Fibroblast"
    ]
    
    # 3. Reference SMC
    ref_smc <- healthy_data[, 
        healthy_data$tier2_merged == "Smooth Muscle"
    ]
    
    # 4. Reference Fibro
    ref_fibro <- healthy_data[, 
        healthy_data$tier2_merged == "Fibroblast"
    ]
    
    cat(sprintf("Inflamed Fibro to SMC: %d cells\n", ncol(inflamed_to_smc)))
    cat(sprintf("Inflamed Fibro to Fibro: %d cells\n", ncol(inflamed_to_fibro)))
    cat(sprintf("Reference SMC: %d cells\n", ncol(ref_smc)))
    cat(sprintf("Reference Fibro: %d cells\n", ncol(ref_fibro)))
    
    # _____________________________________________
    # Find SMC vs Fibroblast markers in Reference
    # _____________________________________________
    
    cat("\n--- Finding SMC vs Fibroblast markers in healthy reference ---\n")
    
    ref_comparison <- cbind(ref_smc, ref_fibro)
    ref_comparison$celltype <- c(
        rep("SMC", ncol(ref_smc)),
        rep("Fibroblast", ncol(ref_fibro))
    )
    
    cat("Computing marker genes...\n")
    ref_markers <- findMarkers(
        x = ref_comparison,
        groups = ref_comparison$celltype,
        pval.type = "all",
        direction = "up"
    )
    
    # Get SMC markers
    smc_markers <- ref_markers[["SMC"]]
    smc_markers <- smc_markers[smc_markers$FDR < 0.01, ]
    all_smc_genes <- rownames(smc_markers)
    
    # Get Fibroblast markers
    fib_markers <- ref_markers[["Fibroblast"]]
    fib_markers <- fib_markers[fib_markers$FDR < 0.01, ]
    all_fib_genes <- rownames(fib_markers)
    
    cat(sprintf("✓ SMC markers found: %d\n", length(all_smc_genes)))
    cat(sprintf("✓ Fibroblast markers found: %d\n", length(all_fib_genes)))
    
    # _________________________
    # Make colData compatible
    # _________________________
    
    cat("\nMaking colData compatible across datasets...\n")
    
    common_cols <- intersect(colnames(colData(inflamed_to_smc)), colnames(colData(ref_smc)))
    
    colData(inflamed_to_smc) <- colData(inflamed_to_smc)[, common_cols]
    colData(ref_smc) <- colData(ref_smc)[, common_cols]
    colData(inflamed_to_fibro) <- colData(inflamed_to_fibro)[, common_cols]
    colData(ref_fibro) <- colData(ref_fibro)[, common_cols]
    
    # ______________________________________
    # Create pseudobulk profiles
    # ______________________________________
    
    cat("\n--- Preparing pseudobulk data ---\n")
    
    # Rename cells for uniqueness
    colnames(inflamed_to_smc) <- paste0("Inflamed_SMC_", seq_len(ncol(inflamed_to_smc)))
    colnames(inflamed_to_fibro) <- paste0("Inflamed_Fibro_", seq_len(ncol(inflamed_to_fibro)))
    colnames(ref_smc) <- paste0("RefSMC_", seq_len(ncol(ref_smc)))
    colnames(ref_fibro) <- paste0("RefFibro_", seq_len(ncol(ref_fibro)))
    
    # Combine all groups
    combined_all <- cbind(inflamed_to_smc, inflamed_to_fibro, ref_smc, ref_fibro)
    
    combined_all$group <- c(
        rep("Inflamed Fibro to SMC", ncol(inflamed_to_smc)),
        rep("Inflamed Fibro to Fibro", ncol(inflamed_to_fibro)),
        rep("Reference SMC", ncol(ref_smc)),
        rep("Reference Fibro", ncol(ref_fibro))
    )
    
    # Create pseudobulk aggregation for all reference marker genes
    cat("Creating pseudobulk profiles...\n")
    
    groups <- unique(combined_all$group)
    all_genes <- unique(c(all_smc_genes, all_fib_genes))
    
    pb_matrix_all <- matrix(0, nrow = length(all_genes), ncol = length(groups))
    rownames(pb_matrix_all) <- all_genes
    colnames(pb_matrix_all) <- groups
    
    for (i in seq_along(groups)) {
        group_cells <- combined_all[, combined_all$group == groups[i]]
        pb_matrix_all[, i] <- rowMeans(assay(group_cells[all_genes, ], "logcounts"))
    }
    
    cat(sprintf("Pseudobulk matrix dimensions: %d genes x %d groups\n", 
                nrow(pb_matrix_all), ncol(pb_matrix_all)))
    
    # _______________________________
    # Find top discriminating genes
    # _______________________________
    
    cat("\n--- Finding top genes separating Inflamed Fibro to SMC from Inflamed Fibro to Fibro ---\n")
    
    expr_diff <- pb_matrix_all[, "Inflamed Fibro to SMC"] - pb_matrix_all[, "Inflamed Fibro to Fibro"]
    top_diff_genes <- names(sort(abs(expr_diff), decreasing = TRUE)[1:50])
    
    # Separate by marker type and get top 10 from each
    top_smc_from_diff <- top_diff_genes[top_diff_genes %in% all_smc_genes]
    top_fib_from_diff <- top_diff_genes[top_diff_genes %in% all_fib_genes]
    
    selected_genes_smc <- head(top_smc_from_diff[!is.na(top_smc_from_diff)], 10)
    selected_genes_fib <- head(top_fib_from_diff[!is.na(top_fib_from_diff)], 10)
    
    cat(sprintf("Top 10 SMC markers: %s\n", paste(selected_genes_smc, collapse = ", ")))
    cat(sprintf("Top 10 Fibro markers: %s\n", paste(selected_genes_fib, collapse = ", ")))
    
    # Split into two matrices by gene type
    pb_smc <- pb_matrix_all[selected_genes_smc, ]
    pb_fibro <- pb_matrix_all[selected_genes_fib, ]
    
    # _______________________
    # Create ComplexHeatmap
    # _______________________
    
    cat("\n--- Creating ComplexHeatmap ---\n")
    
    # Create column annotation (same for both heatmaps)
    col_colors <- c(
        "Inflamed Fibro to SMC" = "#A23B72",
        "Inflamed Fibro to Fibro" = "#2E86AB",
        "Reference SMC" = "#C73E1D",
        "Reference Fibro" = "#F18F01"
    )
    
    col_anno <- HeatmapAnnotation(
        `Cell Group` = factor(colnames(pb_matrix_all), levels = colnames(pb_matrix_all)),
        col = list(`Cell Group` = col_colors),
        annotation_height = unit(0.6, "cm"),
        annotation_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9)
        )
    )
    
    # Determine color scale limits (white to red, one-directional)
    color_min <- min(pb_matrix_all, na.rm = TRUE)
    color_max <- max(pb_matrix_all, na.rm = TRUE)
    
    # Create SMC marker heatmap
    ht_smc <- Heatmap(
        pb_smc,
        name = "Mean\nlogcounts",
        column_title = method_name,
        row_title = "SMC Markers",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_names_side = "bottom",
        row_names_side = "left",
        top_annotation = col_anno,
        cell_fun = function(j, i, x, y, w, h, fill) {
            grid::grid.text(sprintf("%.2f", pb_smc[i, j]), x, y, 
                           gp = grid::gpar(fontsize = 8, col = "black", fontface = "bold"))
        },
        col = colorRamp2(
            c(color_min, color_max),
            c("white", "#b2182b")
        ),
        height = unit(4, "cm"),
        width = unit(5, "cm"),
        heatmap_legend_param = list(
            title = "Mean\nlogcounts",
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9),
            direction = "v",
            legend_height = unit(2, "cm")
        ),
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.5)
    )
    
    # Create Fibroblast marker heatmap
    ht_fibro <- Heatmap(
        pb_fibro,
        name = "Mean\nlogcounts",
        row_title = "Fibroblast Markers",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_names_side = "bottom",
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
            grid::grid.text(sprintf("%.2f", pb_fibro[i, j]), x, y, 
                           gp = grid::gpar(fontsize = 8, col = "black", fontface = "bold"))
        },
        col = colorRamp2(
            c(color_min, color_max),
            c("white", "#b2182b")
        ),
        height = unit(4, "cm"),
        width = unit(5, "cm"),
        heatmap_legend_param = list(
            title = "Mean\nlogcounts",
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 9),
            direction = "v",
            legend_height = unit(2, "cm")
        ),
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.5)
    )
    
    # Combine heatmaps vertically
    ht_combined <- ht_smc %v% ht_fibro
    
    return(ht_combined)
}

# _____________________________________
# Process All Four Annotation Methods
# _____________________________________

cat("\n=== CREATING HEATMAPS FOR ALL METHODS ===\n")

ht_azimuth <- create_annotation_heatmaps("Azimuth", "azimuth_celltype_l1_merged")
ht_singler <- create_annotation_heatmaps("SingleR", "singler_annotations_merged")
ht_celltypist <- create_annotation_heatmaps("CellTypist", "celltypist_predicted_labels_merged")
ht_scarches <- create_annotation_heatmaps("scArches", "scvi_prediction_merged")

# _______________________
# Save Individual Method 
# _______________________

cat("\n--- Saving individual heatmaps ---\n")

# Save Azimuth as main Figure S6
cat("Saving Azimuth heatmap as Fig_S6...\n")
dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)
png("figures/supp/merfish/Fig_S6_Annotation_Patterns_Heatmap_Azimuth.png", 
    width = 5.5, height = 4, res = 300, units = "in")
draw(ht_azimuth, padding = unit(c(1, 1, 1, 1), "mm"))
dev.off()

cat("✓ Azimuth heatmap saved as main Fig_S6\n")