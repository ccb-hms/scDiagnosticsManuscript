# -----------------------------------------------
# MERFISH Mouse Colon IBD - Figure 3 Cell-Cell Interaction (with SingleCellSignalR)
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(SingleCellSignalR)
library(circlize) # For the chord diagram
library(ggplot2)
library(patchwork)

# -----------------------------------------------

# Load the processed DSS9 data object (which should have 'anomaly_status')
dss9_data <- readRDS("data/merfish/dss9_data.rds")
# Re-run anomaly detection from previous scripts if 'anomaly_status' is missing.

# __________________________________________________________________
# 1. Prepare Data for SingleCellSignalR
# __________________________________________________________________
# The setup is very similar to the other tools.

# --- Create the "Interaction Group" Label ---
dss9_data$interaction_group <- as.character(dss9_data$tier2)
dss9_data$interaction_group[grepl("Fibro|IAF", dss9_data$tier2) & dss9_data$anomaly_status == "Anomalous"] <- "Fibroblast_Anomalous"
dss9_data$interaction_group[grepl("Fibro", dss9_data$tier2) & dss9_data$anomaly_status == "Typical"] <- "Fibroblast_Typical"
dss9_data$interaction_group[grepl("Macrophage", dss9_data$tier2)] <- "Macrophage"
dss9_data$interaction_group[grepl("Monocyte", dss9_data$tier2)] <- "Monocyte"
dss9_data$interaction_group[grepl("Neutrophil", dss9_data$tier2)] <- "Neutrophil"

main_groups <- c("Fibroblast_Anomalous", "Fibroblast_Typical", "Macrophage", "Monocyte", "Neutrophil")
dss9_data$interaction_group[!dss9_data$interaction_group %in% main_groups] <- "Other"
dss9_data$interaction_group <- factor(dss9_data$interaction_group)

# --- Prepare the inputs for the main function ---
# SingleCellSignalR requires a raw counts matrix and a vector of cell labels.
counts_matrix <- as.matrix(counts(dss9_data))
cell_labels <- dss9_data$interaction_group
names(cell_labels) <- colnames(dss9_data)


# __________________________________________________________________
# 2. Run Cell-Cell Interaction Inference with SingleCellSignalR
# __________________________________________________________________
cat("\nRunning SingleCellSignalR to infer cell-cell interactions...\n")

# This is the main function call. It performs all calculations internally.
# We specify 'mouse' as the species.
signal <- cell_signaling(
    data = counts_matrix,
    genes = rownames(counts_matrix),
    cluster = cell_labels,
    species = "mouse"
)

cat("SingleCellSignalR analysis complete.\n")


# __________________________________________________________________
# 3. Create the Visualization (A Two-Panel Figure)
# __________________________________________________________________
# We will use SingleCellSignalR's built-in plotting functions.

# --- Panel A: The Overall Interaction Score Heatmap ---
# This is SingleCellSignalR's specialty. It shows a heatmap of the overall 
# communication "score" from each sender to each receiver.
# We will use this to show that 'Fibroblast_Anomalous' is a major signaling hub.

# The function call is simple, but we will save the plot object to customize it.
panel_a_plot <- signal_heatmap(
    signal, 
    font.size = 10
)

# Customize the ggplot object
panel_a <- panel_a_plot + 
    labs(title = "A) Overall Interaction Strength",
         subtitle = "'Fibroblast_Anomalous' is a major signaling hub") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))


# --- Panel B: A "Spotlight" Chord Diagram on a Key Ligand-Receptor Pair ---
# Let's focus on the CXCL-CXCR interaction, which is critical for neutrophil recruitment.
# We will find a specific pair, e.g., Cxcl1 -> Cxcr2.

# First, extract the data for the chord diagram
# The 'signal$links' data frame contains all the inferred interactions.
lr_links <- signal$links
# Filter for our pathway of interest
cxcl_links <- lr_links[grepl("Cxcl", lr_links$ligand) & grepl("Cxcr", lr_links$receptor), ]

# We need to manually create the plot using the 'circlize' package, as SingleCellSignalR's
# built-in chord diagram plots all interactions, which is too messy.
# Prepare the data for the chord diagram: a matrix of interaction scores.
chord_data <- cxcl_links %>%
    filter(receptor.cluster == "Neutrophil",
           ligand.cluster %in% c("Fibroblast_Anomalous", "Fibroblast_Typical")) %>%
    group_by(ligand.cluster, receptor.cluster) %>%
    summarise(score = sum(score), .groups = 'drop')

# Create the chord diagram
# We have to capture the plot output using recordPlot() because circlize doesn't return a ggplot object.
pdf(NULL) # Prevent circlize from opening a new graphics device
circos.clear()
chordDiagram(chord_data, transparency = 0.5)
title("B) Spotlight: CXCL -> CXCR Signaling to Neutrophils")
panel_b_recorded <- recordPlot()
dev.off()


# --- Combine and Save the Final Figure ---
# We use patchwork to combine the ggplot (Panel A) and the recorded base R plot (Panel B).
final_interaction_figure <- panel_a + panel_b_recorded +
    plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
        title = "Figure 3: Anomalous fibroblasts adopt a pro-inflammatory communication profile",
        theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )

print(final_interaction_figure)

# Save the figure
ggsave("MERFISH_Figure3_Interactions_SignalR.pdf", 
       plot = final_interaction_figure, 
       width = 14, 
       height = 7, 
       units = "in")