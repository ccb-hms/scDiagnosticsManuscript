# ----------------------------------------------------
# COVID-19 PBMC - UMAP Visualization (scVI/scArches)
# ----------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(pals) 
library(viridis)
library(scattermore) 

# ----------------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("COVID-19 UMAP VISUALIZATION (scVI Results)\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading processed and annotated datasets...\n")
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# Verify that the required UMAP coordinates exist
if (!"UMAP_scVI" %in% reducedDimNames(normal_data) || !"UMAP_scVI" %in% reducedDimNames(covid_data)) {
    stop("Error: 'UMAP_scVI' not found. Please re-run the scVI pipeline.")
}
cat("✓ 'UMAP_scVI' coordinates found in both datasets.\n")

# _____________________________________________
# Create a Single, Combined Plotting Data Frame
# _____________________________________________

cat("\nCombining all data into a single data frame for plotting...\n")

# Extract UMAP coordinates
umap_normal <- as.data.frame(reducedDim(normal_data, "UMAP_scVI"))
colnames(umap_normal) <- c("UMAP1", "UMAP2")
umap_covid <- as.data.frame(reducedDim(covid_data, "UMAP_scVI"))
colnames(umap_covid) <- c("UMAP1", "UMAP2")

# Extract and prepare metadata
meta_normal <- as.data.frame(colData(normal_data)) %>%
    select(author_cell_type_merged) %>%
    mutate(disease_status = "Healthy",
           hybrid_cell_type = author_cell_type_merged)

meta_covid <- as.data.frame(colData(covid_data)) %>%
    select(author_cell_type_merged, azimuth_celltype_l1_merged) %>%
    mutate(disease_status = "COVID-19",
           hybrid_cell_type = azimuth_celltype_l1_merged)

# Combine into one data frame
plot_df <- bind_rows(
    bind_cols(umap_normal, meta_normal),
    bind_cols(umap_covid, meta_covid)
)
cat(sprintf("✓ Combined plotting data frame created with %d cells.\n", nrow(plot_df)))

# __________________________
# IFN Signature Computation
# __________________________

cat("\nComputing and adding IFN signature scores...\n")
common_genes <- intersect(rownames(normal_data), rownames(covid_data))
combined_logcounts <- cbind(
    assay(normal_data, "logcounts")[common_genes, ],
    assay(covid_data, "logcounts")[common_genes, ]
)

ifn_genes <- c("BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", "IFI6", "IFIT3", 
               "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", "PARP9", "PLSCR1", "SAMD9", 
               "SAMD9L", "SP110", "STAT1", "TRIM22", "UBE2L6", "XAF1", "IRF7")
available_ifn_genes <- intersect(ifn_genes, rownames(combined_logcounts))

if(length(available_ifn_genes) > 0) {
    ifn_expression <- combined_logcounts[available_ifn_genes, , drop = FALSE]
    plot_df$IFN_signature <- colMeans(ifn_expression)
    cat("✓ IFN signature score computed and added to plotting data frame.\n")
} else {
    plot_df$IFN_signature <- 0
}

# ___________________
# UMAP Visualization
# ___________________

cat("\nCreating UMAP visualizations...\n")

cell_type_colors <- c(
    
    # T Cell & ILC Lineage
    "CD4 T"       = "#8A2BE2",  
    "CD8 T"       = "#FF8C00",  
    "Treg"        = "#DA70D6",  
    "NKT"         = "#9932CC",  
    "MAIT"        = "#FF4500",  
    "gdT"         = "#FFA500",  
    "ILC"         = "#FF69B4",  
    
    # B Cell & Plasma Lineage
    "B cell"      = "#20B2AA",  
    "Plasmablast" = "#32CD32",  
    "Plasma cell" = "#006400",  

    # Myeloid & Monocyte Lineage
    "CD14 mono"   = "#008080",  
    "CD16 mono"   = "#1E90FF",  
    "Mono_prolif" = "#ADD8E6",  
    "DC1"         = "#A0522D",  
    "DC2"         = "#D2691E",  
    "DC3"         = "#F4A460",  
    "ASDC"        = "#8B4513",  
    "pDC"         = "#CD853F",  
    "DC_prolif"   = "#FFE4B5",  
    
    # NK Cell Lineage
    "NK_16hi"     = "#B22222",  
    "NK_56hi"     = "#8B4513",  
    "NK_prolif"   = "#FFC0CB",  
    
    # Other Types
    "HSPC"        = "#F0E68C",  
    "Platelets"   = "#BDB76B",  
    "RBC"         = "#E9967A"   
)

# Define the disease status palette
disease_colors <- c("Healthy" = "#5DADE2", "COVID-19" = "#FFB6C1") 

cat("✓ Custom color palettes created.\n")

# Create a Reusable ggplot Theme
publication_theme <- theme_classic(base_size = 12) +
    theme(
        plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.line = element_line(colour = "black", linewidth = 0.5)
    )

# Define point aesthetics for consistency
point_size <- 0.1
point_alpha <- 0.3

# Plot 1: UMAP colored by disease status
cat("Creating UMAP plot by disease status...\n")  
p_disease <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = disease_status)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = disease_colors) + 
    labs(color = "Condition", x = "UMAP 1", y = "UMAP 2") +
    publication_theme +
    guides(color = guide_legend(override.aes = list(size = 5)))

# Plot 2: UMAP colored by hybrid cell types
cat("Creating UMAP plot by hybrid cell type...\n")
label_df <- plot_df %>%
    filter(!is.na(hybrid_cell_type)) %>%
    group_by(hybrid_cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = 'drop')

p_celltype <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = hybrid_cell_type), size = point_size, alpha = point_alpha + 0.1) +
    scale_color_manual(values = cell_type_colors, name = "Cell Type", na.value = "grey80") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) 

# Plot 3: UMAP colored by IFN signature
cat("Creating UMAP plot by IFN signature...\n")
p_ifn <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = IFN_signature)) +
    geom_point(size = point_size, alpha = 0.5) +
    scale_color_viridis_c(name = "IFN Score", option = "plasma") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    publication_theme

cat("✓ UMAP plots created.\n")

# _______________________
# Display and Save Plots
# _______________________

cat("\nDisplaying and saving plots in order: Disease, Cell Type, IFN Signature...\n")
print(p_disease)
print(p_celltype)
print(p_ifn)

dir.create("figures/covid", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/covid/scvi_umap_1_by_disease.png", p_disease, width = 10, height = 8, dpi = 600)
ggsave("figures/covid/scvi_umap_2_by_hybrid_celltype.png", p_celltype, width = 10, height = 8, dpi = 600)
ggsave("figures/covid/scvi_umap_3_by_ifn.png", p_ifn, width = 10, height = 8, dpi = 600)

cat("✓ Plots saved to figures/covid/ directory.\n")
cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("scVI UMAP VISUALIZATION COMPLETED\n")
cat(paste(rep("=", 45), collapse = ""), "\n")
