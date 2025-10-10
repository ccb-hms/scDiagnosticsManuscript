# -----------------------------------------------
# MERFISH Mouse Colon IBD - Figure 3 Cell-Cell Interaction (with Liana)
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(liana)
library(ggplot2)
library(patchwork)

# -----------------------------------------------

# Load the processed DSS9 data object (which should have 'anomaly_status')
dss9_data <- readRDS("data/merfish/dss9_data.rds")
# Re-run anomaly detection from previous scripts if 'anomaly_status' is missing.

# __________________________________________________________________
# 1. Prepare Data for Liana
# __________________________________________________________________
# This part is very similar to the CellChat prep.

# --- Create the "Interaction Group" Label ---
dss9_data$interaction_group <- as.character(dss9_data$tier2)
dss9_data$interaction_group[grepl("Fibro|IAF", dss9_data$tier2) & dss9_data$anomaly_status == "Anomalous"] <- "Fibroblast (Anomalous)"
dss9_data$interaction_group[grepl("Fibro", dss9_data$tier2) & dss9_data$anomaly_status == "Typical"] <- "Fibroblast (Typical)"
dss9_data$interaction_group[grepl("Macrophage", dss9_data$tier2)] <- "Macrophage"
dss9_data$interaction_group[grepl("Monocyte", dss9_data$tier2)] <- "Monocyte"
dss9_data$interaction_group[grepl("Neutrophil", dss9_data$tier2)] <- "Neutrophil"

# Simplify and group remaining types
main_groups <- c("Fibroblast (Anomalous)", "Fibroblast (Typical)", "Macrophage", "Monocyte", "Neutrophil")
dss9_data$interaction_group[!dss9_data$interaction_group %in% main_groups] <- "Other"
dss9_data$interaction_group <- factor(dss9_data$interaction_group)

cat("Created 'interaction_group' labels. Distribution:\n")
print(table(dss9_data$interaction_group))


# __________________________________________________________________
# 2. Run Cell-Cell Interaction Inference with Liana
# __________________________________________________________________

# Liana works directly on the SingleCellExperiment object.
# We will run the 'Connectome' method, which is a powerful consensus-based approach.
cat("\nRunning Liana to infer cell-cell interactions...\n")
liana_results <- liana_wrap(
    dss9_data,
    assay_name = "logcounts", # Liana can work with logcounts
    method = "connectome",   # A robust consensus method
    resource = "mouse",    # Specify the mouse ligand-receptor database
    group_by = "interaction_group"
)

# The result is a data frame of all predicted interactions.
cat("Liana analysis complete. Found", nrow(liana_results), "potential interactions.\n")


# __________________________________________________________________
# 3. Create the Visualization (A Two-Panel Figure)
# __________________________________________________________________
# We will manually create a differential dot plot from the liana_results data frame.

# --- Step 3a: Perform Differential Analysis ---
# We want to see which interactions are specific to our 'Anomalous' vs 'Typical' senders.
# We filter for interactions originating from our two fibroblast groups.
fibro_interactions <- liana_results %>%
    filter(source %in% c("Fibroblast (Anomalous)", "Fibroblast (Typical)"))

# Let's focus on the key inflammatory ligands from the paper.
key_ligands <- c("Il11", "Il1b", "Mmp3", "Mmp10", "Cxcl1", "Cxcl5", "Csf3")

# Filter for interactions involving these key ligands and the most important targets
differential_interactions <- fibro_interactions %>%
    filter(ligand.expr_req %in% key_ligands,
           target %in% c("Neutrophil", "Macrophage", "Monocyte", "Fibroblast (Anomalous)", "Fibroblast (Typical)")) %>%
    # Use a ranking score to find the most important interactions
    # We use 'magnitude_rank' which combines expression and specificity
    arrange(magnitude_rank) %>%
    # Select the top interactions for each sender group to keep the plot clean
    group_by(source) %>%
    top_n(15, wt = -magnitude_rank) # Negative sign because lower rank is better

cat("Filtered down to", nrow(differential_interactions), "key differential interactions for plotting.\n")


# --- Step 3b: Create the Differential Dot Plot (Panel A) ---
# This plot will show the strength of key inflammatory signals from each fibroblast group.
panel_a <- ggplot(differential_interactions, 
                  aes(x = target, y = ligand.expr_req, size = specificity_rank, color = aggregate_rank)) +
    geom_point() +
    facet_wrap(~source, ncol = 2) + # This creates the side-by-side comparison
    scale_color_viridis_c(direction = -1, name = "Aggregate Rank") +
    scale_size(range = c(6, 1), name = "Specificity Rank") + # Invert size scale so better ranks are bigger
    labs(
        title = "A) Gained Pro-inflammatory Signaling from Anomalous Fibroblasts",
        subtitle = "Comparison of outgoing signals from 'Typical' vs 'Anomalous' fibroblasts",
        x = "Target Cell Type",
        y = "Ligand Gene"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold")
    )


# --- Step 3c: Create a "Spotlight" Chord Diagram (Panel B) ---
# Liana provides a convenient function for this. We will spotlight the CXCL signaling.
# First, filter the results for just the CXCL pathway interactions.
cxcl_interactions <- liana_results %>% 
    filter(grepl("CXCL", interaction))

panel_b <- liana_chord(
    cxcl_interactions, 
    source_groups = c("Fibroblast (Anomalous)", "Fibroblast (Typical)"),
    target_groups = c("Neutrophil", "Macrophage")
) + 
    labs(title = "B) Spotlight on CXCL Pathway", subtitle = "Signals to Neutrophils and Macrophages") +
    theme(plot.title = element_text(hjust = 0.5, face="bold"), 
          plot.subtitle = element_text(hjust = 0.5))


# --- Combine and Save the Final Figure ---
final_interaction_figure <- panel_a / panel_b + 
    plot_layout(heights = c(2, 1.5)) +
    plot_annotation(
        title = "Figure 3: Anomalous fibroblasts adopt a pro-inflammatory communication profile",
        theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )

print(final_interaction_figure)

# Save the figure
ggsave("MERFISH_Figure3_Interactions_Liana.pdf", 
       plot = final_interaction_figure, 
       width = 12, 
       height = 10, 
       units = "in")