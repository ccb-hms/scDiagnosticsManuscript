# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 2 Code
# -----------------------------------------------

library(SingleCellExperiment)
library(ggplot2)

# -----------------------------------------------

# __________
# Load data
# __________

covid_data <- readRDS("data/covid/covid_data_sce.rds")
normal_data <- readRDS("data/covid/normal_data_sce.rds")

# _______________________
# Calculate IFN Signature
# _______________________

yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

# Find available genes
common_genes <- intersect(rownames(covid_data), rownames(normal_data))
signature_genes_avail <- intersect(yoshida_ifn_signature, common_genes)

# Calculate IFN score as mean expression of signature genes
ifn_expr <- assay(covid_data, "logcounts")[signature_genes_avail, ]
ifn_score <- colMeans(ifn_expr, na.rm = TRUE)

# Extract UMAP coordinates
umap_coords <- reducedDim(covid_data, "UMAP_scVI")

# ____________
# Setup Theme
# ____________

theme_set(theme_minimal() + theme(
  axis.title = element_text(size = 11, face = "bold"),
  axis.text = element_text(size = 10),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 9, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 10, face = "bold")
))

# ______________________________
# Figure S2: IFN Signature UMAP
# ______________________________

umap_data_ifn <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    IFN_Score = ifn_score
)

fig_s2 <- ggplot(umap_data_ifn, aes(x = UMAP1, y = UMAP2, color = IFN_Score)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_viridis_c(option = "B", name = "IFN Score") +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("IFN Signature Score (Query Cells)") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 11, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 10, color = "#4B5563", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold", color = "#374151"),
        legend.text = element_text(size = 9, color = "#4B5563"),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.5, fill = NA),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA),
        aspect.ratio = 1
    )

ggsave("figures/supp/covid/Fig_S2_ifn_signature.png", fig_s2, width = 8, height = 8, dpi = 600)

print("Supplementary Figure S2 saved")

# ________
# Summary
# ________

print("Supplementary Figure S2 (IFN Signature) complete!")
print("Saved in: figures/supp/covid/")