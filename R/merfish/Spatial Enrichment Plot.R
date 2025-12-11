# --------------------------------------------------
# MERFISH Mouse Colon IBD - Spatial Enrichment Figure
# --------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scDiagnostics)
library(FNN) 

# --------------------------------------------------

# __________
# Load Data
# __________

cat("--- Loading data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Setting the seed
set.seed(0)

# __________________
# Anomaly Detection
# __________________

cat("--- Running Anomaly Detection ---\n")
is_fibro_lineage <- function(x) { grepl("^Fibro|^IAF", x) }
healthy_data$analysis_class <- ifelse(is_fibro_lineage(healthy_data$tier2), "Fibroblast_Lineage", "Other")
dss9_data$analysis_class <- ifelse(is_fibro_lineage(dss9_data$tier2), "Fibroblast_Lineage", "Other")

anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "analysis_class", 
    ref_cell_type_col = "analysis_class", 
    cell_types = "Fibroblast_Lineage",
    pc_subset = 1:4, 
    anomaly_threshold = 0.5,
    max_cells_query = NULL,
    max_cells_ref = NULL
)

anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# ______________________
# Spatial Calculations
# ______________________

spatial_df <- data.frame(
    barcode = colnames(dss9_data),
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    broad_class = dss9_data$analysis_class
)

# Define Targets
coords_neutro   <- spatial_df |> filter(grepl("Neutrophil", ground_truth_type, ignore.case=TRUE)) |> select(x,y)
coords_inflamed <- spatial_df |> filter(grepl("IAE|Epithelial \\(Clu\\+\\)", ground_truth_type, ignore.case=TRUE)) |> select(x,y)
coords_stem     <- spatial_df |> filter(grepl("Stem|TA", ground_truth_type, ignore.case=TRUE)) |> select(x,y)

# Calculate Distances for Fibroblasts Only
fibro_df <- spatial_df |> filter(broad_class == "Fibroblast_Lineage")
fibro_df$prediction <- ifelse(fibro_df$barcode %in% anomalous_barcodes, "Anomalous", "Typical")

get_dists <- function(coords) {
    if(nrow(coords) == 0) return(rep(NA, nrow(fibro_df)))
    FNN::get.knnx(data = as.matrix(coords), query = as.matrix(fibro_df[, c("x", "y")]), k = 1)$nn.dist
}

fibro_df$dist_neutro   <- get_dists(coords_neutro)
fibro_df$dist_inflamed <- get_dists(coords_inflamed)
fibro_df$dist_stem     <- get_dists(coords_stem)

# _________________________________
# Enrichment Calculation (Binning)
# _________________________________

# 1. Calculate Statistics for Title (Wilcoxon)
calc_pval_str <- function(col_name) {
    p <- wilcox.test(fibro_df[[col_name]] ~ fibro_df$prediction)$p.value
    if(p < 2.2e-16) return("p < 2.2e-16")
    return(paste0("p = ", format.pval(p, digits=2)))
}

# Define Panel Titles
lbl_neutro   <- paste0("Neutrophils (Acute)\n", calc_pval_str("dist_neutro"))
lbl_inflamed <- paste0("Inflamed Epithelium (Damage)\n", calc_pval_str("dist_inflamed"))
lbl_stem     <- paste0("Stem/Crypt (Homeostasis)\n", calc_pval_str("dist_stem"))

# 2. Reshape to Long Format
long_df <- fibro_df |>
    select(barcode, prediction, dist_neutro, dist_inflamed, dist_stem) |>
    tidyr::pivot_longer(
        cols = c(dist_neutro, dist_inflamed, dist_stem),
        names_to = "target_type",
        values_to = "distance"
    ) |>
    mutate(panel_label = case_when(
        target_type == "dist_neutro"   ~ lbl_neutro,
        target_type == "dist_inflamed" ~ lbl_inflamed,
        target_type == "dist_stem"     ~ lbl_stem
    ))

# 3. Create Distance Bins (Focus on 0-100um)
long_df$dist_bin <- cut(long_df$distance, 
                        breaks = c(0, 25, 50, 75, 100, 150, Inf),
                        labels = c("0-25", "25-50", "50-75", "75-100", "100-150", ">150"))

# 4. Calculate Enrichment Ratios
enrichment_df <- long_df |>
    # Count cells per bin per type
    group_by(panel_label, dist_bin, prediction) |>
    summarise(count = n(), .groups = "drop") |>
    
    # Normalize by total library size (Typical vs Anomalous)
    group_by(prediction) |>
    mutate(total_cells = sum(count)) |>
    ungroup() |>
    mutate(prop = count / total_cells) |>
    
    # Pivot to get Anomalous vs Typical side-by-side
    select(-count, -total_cells) |>
    tidyr::pivot_wider(names_from = prediction, values_from = prop) |>
    
    # Calculate Log2 Fold Change
    mutate(log2fc = log2(Anomalous / Typical)) |>
    
    # Filter out the ">150" bin as it's just background noise
    filter(dist_bin != ">150")

# _________________________________________
# *** FORCE PANEL ORDER HERE ***
# _________________________________________

enrichment_df$panel_label <- factor(enrichment_df$panel_label, 
                                    levels = c(lbl_neutro, lbl_inflamed, lbl_stem))

# __________________
# Generate Plot
# __________________

# Custom Color Logic
enrichment_df$enrichment_class <- ifelse(enrichment_df$log2fc > 0, "Anomalous Enriched", "Typical Enriched")

p <- ggplot(enrichment_df, aes(x = dist_bin, y = log2fc, fill = enrichment_class)) +
    geom_col(width = 0.7, color = "black", size = 0.3) +
    
    # Horizontal line at 0 (No difference)
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    
    # Facets
    facet_wrap(~panel_label, scales = "fixed") +
    
    # Colors
    scale_fill_manual(values = c("Anomalous Enriched" = "#D9230F", "Typical Enriched" = "#228B22")) +
    
    # Labels
    labs(
        title = "Spatial Enrichment Analysis: The Inflammatory Niche",
        subtitle = "Log2 Fold Enrichment of Anomalous vs. Typical Fibroblasts at specific distances",
        x = expression(bold("Distance from Feature ("*mu*"m)")),
        y = expression(bold("Log"[2]*" Enrichment (Anomalous / Typical)")),
        fill = "Enrichment Status"
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
    )

print(p)

ggsave("figures/merfish/spatial_enrichment_analysis.png", plot = p, width = 12, height = 7, dpi = 600)
cat("✓ Spatial Enrichment Bar Plot saved.\n")