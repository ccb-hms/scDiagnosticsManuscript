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
    
    inflamed_to_smc <- dss9_data[, dss9_data$tier2_merged == "Inflamed Fibroblast" & dss9_data[[query_annot_col]] == "Smooth Muscle"]
    inflamed_to_fibro <- dss9_data[, dss9_data$tier2_merged == "Inflamed Fibroblast" & dss9_data[[query_annot_col]] == "Fibroblast"]
    ref_smc <- healthy_data[, healthy_data$tier2_merged == "Smooth Muscle"]
    ref_fibro <- healthy_data[, healthy_data$tier2_merged == "Fibroblast"]
    
    # _____________________________________________
    # Find SMC vs Fibroblast markers in Reference
    # _____________________________________________
    
    ref_comparison <- cbind(ref_smc, ref_fibro)
    ref_comparison$celltype <- c(rep("SMC", ncol(ref_smc)), rep("Fibroblast", ncol(ref_fibro)))
    
    ref_markers <- findMarkers(x = ref_comparison, groups = ref_comparison$celltype, pval.type = "all", direction = "up")
    
    smc_markers <- ref_markers[["SMC"]]
    all_smc_genes <- rownames(smc_markers[smc_markers$FDR < 0.01, ])
    
    fib_markers <- ref_markers[["Fibroblast"]]
    all_fib_genes <- rownames(fib_markers[fib_markers$FDR < 0.01, ])
    
    # _________________________
    # Make colData compatible
    # _________________________
    
    common_cols <- intersect(colnames(colData(inflamed_to_smc)), colnames(colData(ref_smc)))
    colData(inflamed_to_smc) <- colData(inflamed_to_smc)[, common_cols]
    colData(ref_smc) <- colData(ref_smc)[, common_cols]
    colData(inflamed_to_fibro) <- colData(inflamed_to_fibro)[, common_cols]
    colData(ref_fibro) <- colData(ref_fibro)[, common_cols]
    
    # ______________________________________
    # Create pseudobulk profiles
    # ______________________________________
    
    # Cleaner, more concise column names
    group_names <- c(
        "Misannotated\n(as SMC)",
        "Reference\n(True SMC)",
        "Annotated\n(as Fibroblast)",
        "Reference\n(True Fibroblast)"
    )
    
    combined_all <- cbind(inflamed_to_smc, ref_smc, inflamed_to_fibro, ref_fibro)
    combined_all$group <- factor(
        c(
            rep(group_names[1], ncol(inflamed_to_smc)),
            rep(group_names[2], ncol(ref_smc)),
            rep(group_names[3], ncol(inflamed_to_fibro)),
            rep(group_names[4], ncol(ref_fibro))
        ), 
        levels = group_names
    )
    
    all_genes <- unique(c(all_smc_genes, all_fib_genes))
    pb_matrix_all <- matrix(0, nrow = length(all_genes), ncol = length(group_names))
    rownames(pb_matrix_all) <- all_genes
    colnames(pb_matrix_all) <- group_names
    
    for (i in seq_along(group_names)) {
        group_cells <- combined_all[, combined_all$group == group_names[i]]
        pb_matrix_all[, i] <- rowMeans(assay(group_cells[all_genes, ], "logcounts"))
    }
    
    # _______________________________
    # Find top discriminating genes
    # _______________________________
    
    expr_diff <- pb_matrix_all[, group_names[1]] - pb_matrix_all[, group_names[3]]
    top_diff_genes <- names(sort(abs(expr_diff), decreasing = TRUE)[1:50])
    
    top_smc_from_diff <- top_diff_genes[top_diff_genes %in% all_smc_genes]
    top_fib_from_diff <- top_diff_genes[top_diff_genes %in% all_fib_genes]
    
    selected_genes_smc <- head(top_smc_from_diff[!is.na(top_smc_from_diff)], 10)
    selected_genes_fib <- head(top_fib_from_diff[!is.na(top_fib_from_diff)], 10)
    
    # Raw logcounts matrices
    pb_smc <- pb_matrix_all[selected_genes_smc, ]
    pb_fibro <- pb_matrix_all[selected_genes_fib, ]
    
    # _______________________
    # Create ComplexHeatmap
    # _______________________
    
    cat("\n--- Creating ComplexHeatmap ---\n")
    
    # Establish a clean global color scale
    color_max <- max(c(pb_smc, pb_fibro), na.rm = TRUE)
    # Round up to nearest 0.5 for a clean legend
    color_max <- ceiling(color_max * 2) / 2 
    
    global_col_fun <- colorRamp2(c(0, color_max), c("white", "#b2182b"))
    
    # Text color logic (White text on dark red, Black text on light red/white)
    cell_text_fun <- function(matrix_data) {
        function(j, i, x, y, w, h, fill) {
            val <- matrix_data[i, j]
            text_color <- ifelse(val > (color_max * 0.6), "white", "black")
            grid::grid.text(sprintf("%.2f", val), x, y, 
                           gp = grid::gpar(fontsize = 16, col = text_color, fontface = "bold"))
        }
    }
    
    # SMC Heatmap (Top)
    ht_smc <- Heatmap(
        matrix = pb_smc, 
        name = "Mean\nLogcounts",
        column_title = paste0(method_name, " Annotations on Inflamed Fibroblasts"),
        column_title_gp = gpar(fontsize = 26, fontface = "bold"),
        row_title = "SMC Markers",
        row_title_gp = gpar(fontsize = 22, fontface = "bold"),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE, 
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 18, fontface = "italic"),
        cell_fun = cell_text_fun(pb_smc),
        col = global_col_fun,
        height = unit(10, "cm"),
        width = unit(14, "cm"),
        heatmap_legend_param = list(
            title = "Mean\nLogcounts",
            title_position = "topleft", # Fixes overlap by expanding text rightward
            title_gp = gpar(fontsize = 18, fontface = "bold"),
            labels_gp = gpar(fontsize = 16),
            legend_height = unit(6, "cm"),
            grid_width = unit(1.0, "cm") # Thicker color bar
        ),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1)
    )
    
    # Fibroblast Heatmap (Bottom)
    ht_fibro <- Heatmap(
        matrix = pb_fibro,
        name = "Mean\nLogcounts",
        row_title = "Fibroblast Markers",
        row_title_gp = gpar(fontsize = 22, fontface = "bold"),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = TRUE, 
        column_names_side = "bottom",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 20, fontface = "bold"),
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 18, fontface = "italic"),
        cell_fun = cell_text_fun(pb_fibro),
        col = global_col_fun,
        height = unit(10, "cm"),
        width = unit(14, "cm"),
        show_heatmap_legend = FALSE, # Prevent duplicate legends from clashing
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1)
    )
    
    # Combine vertically
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
    width = 13, height = 12, res = 300, units = "in")

# Added extra padding to the RIGHT (c(bottom, left, top, right)) to give the legend plenty of space
draw(ht_azimuth, padding = unit(c(1, 1, 1, 1.5), "in"), ht_gap = unit(0.5, "cm")) 

dev.off()

cat("✓ Azimuth heatmap saved as main Fig_S6\n")