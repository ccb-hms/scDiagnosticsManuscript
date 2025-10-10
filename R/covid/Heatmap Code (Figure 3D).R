# --------------------------------------
# COVID-19 - Heatmaps Code (Figure 3D)
# --------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)

# --------------------------------------

# Load the processed SCE objects
normal_data <- readRDS("data/covid/normal_data.rds")
covid_data <- readRDS("data/covid/covid_data.rds")

# _________________
# Data Preparation
# _________________

# --- Define the gene groups for plotting ---
genes_nkt <- c("NKG7", "KLRD1", "CD3D", "GZMB")
genes_mono <- c("CD14", "LYZ", "VCAN")
genes_inflam <- c("S100A8", "S100A9", "S100A12", "ISG15", "MX1", "IFITM3")
genes_of_interest <- c(genes_mono, genes_nkt, genes_inflam) # Reordered for better story

# Create a factor for splitting rows. Order matters.
gene_groups <- factor(
  rep(c("Monocyte Markers", "NK T Markers", "Inflammatory Markers"), 
      times = c(length(genes_mono), length(genes_nkt), length(genes_inflam))),
  levels = c("Monocyte Markers", "NK T Markers", "Inflammatory Markers")
)

# Safety check and get indices
genes_present <- genes_of_interest[genes_of_interest %in% rowData(normal_data)$gene_symbol]
gene_indices <- match(genes_present, rowData(normal_data)$gene_symbol)
gene_groups <- gene_groups[genes_of_interest %in% genes_present]

# _________________________________________
# Figure 1 v2.0: Building the 4-Column Pseudo-bulk Matrix
# _________________________________________

# --- Calculate average expression for each of the FOUR groups ---

# Group 1: Healthy Monocytes
idx_mono_healthy <- which(normal_data$cell_type == "CD14-positive monocyte")
pb_mono_healthy <- rowMeans(as.matrix(logcounts(normal_data)[gene_indices, idx_mono_healthy]))

# Group 2: COVID Monocytes
idx_mono_covid <- which(covid_data$cell_type == "CD14-positive monocyte")
pb_mono_covid <- rowMeans(as.matrix(logcounts(covid_data)[gene_indices, idx_mono_covid]))

# Group 3: Healthy NK T cells (The NEW crucial baseline)
idx_nkt_healthy <- which(normal_data$cell_type == "mature NK T cell")
pb_nkt_healthy <- rowMeans(as.matrix(logcounts(normal_data)[gene_indices, idx_nkt_healthy]))

# Group 4: Misannotated NK T cells in COVID
idx_nkt_misannotated <- which(covid_data$cell_type == "mature NK T cell" & 
                              covid_data$singler_annotations == "CD14-positive monocyte")
pb_nkt_misannotated <- rowMeans(as.matrix(logcounts(covid_data)[gene_indices, idx_nkt_misannotated]))


# --- Combine into a single matrix, scale it, and VERIFY ---
pseudo_bulk_matrix <- cbind(
    `Monocyte (Healthy)` = pb_mono_healthy,
    `Monocyte (COVID)` = pb_mono_covid,
    `NK T (Healthy)` = pb_nkt_healthy,
    `Misannotated NK T (COVID)` = pb_nkt_misannotated
)
rownames(pseudo_bulk_matrix) <- genes_present

# !!! IMPORTANT: Look at this matrix to verify the "weird" colors !!!
print("Unscaled Average Log-Counts:")
print(round(pseudo_bulk_matrix, 2))

# Now scale for visualization
scaled_matrix <- t(scale(t(pseudo_bulk_matrix)))

# _________________________________________________________________
# Figure 1 v2.0: Drawing the Improved Heatmap
# _________________________________________________________________

# --- Define the top annotation bar for the four groups ---
# This part remains the same
group_colors <- c("Monocyte (Healthy)" = "steelblue", "Monocyte (COVID)" = "orange",
                  "NK T (Healthy)" = "seagreen", "Misannotated NK T (COVID)" = "firebrick")

top_annotation <- HeatmapAnnotation(
    Group = names(group_colors),
    col = list(Group = group_colors),
    annotation_name_side = "left",
    show_legend = FALSE
)

# --- Draw the final heatmap with spacing adjustments ---
# I've added two new arguments: 'row_gap' and 'padding'

ht_v2_spaced <- Heatmap(
    scaled_matrix, 
    name = "Scaled Expression",
    top_annotation = top_annotation,
    
    # Column settings
    cluster_columns = FALSE,
    column_order = names(group_colors),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 11),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    
    # Row settings
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_split = gene_groups,
    row_names_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 10, fontface = "italic"),
    row_title_rot = 0,
    
    # --- KEY VISUAL CHANGE #1 ---
    row_gap = unit(3, "mm"), # Increase the gap between the gene blocks
    
    # Aesthetics
    col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))
)

# Create the legend object
lgd_group <- Legend(labels = names(group_colors), title = "Group", 
                    legend_gp = gpar(fill = group_colors))

# --- Draw the plot with KEY VISUAL CHANGE #2 ---
# The 'padding' argument adds space around the entire plot.
# The format is: bottom, left, top, right.
# We will add 5mm of padding to the top.
draw(
    ht_v2_spaced, 
    heatmap_legend_list = list(lgd_group), 
    merge_legend = TRUE,
    padding = unit(c(2, 2, 5, 2), "mm") # Add 5mm padding to the top
)
