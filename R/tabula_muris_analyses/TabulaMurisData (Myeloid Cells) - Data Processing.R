# ------------------------------------------
# Process TabulaMurisData - Myeloid Cells
# Option 1: Marrow with Random Split
# ------------------------------------------

# Load libraries
library(TabulaMurisData)
library(SingleCellExperiment)
library(ExperimentHub)
library(scater)
library(scran)
library(SingleR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(knitr)

# source files
source("R/addCellLineageTabulaMuris.R")

# ------------------------------------------

# Set the seed
set.seed(0)

# _____________
# Data loading
# _____________

eh <- ExperimentHub()
tm.facs <- eh[["EH1617"]]  # FACS-sorted Tabula Muris data

# __________________
# Add cell lineage
# __________________

tm.facs <- addCellLineageTabulaMuris(tm.facs)

# Generate count matrix: tissues (rows) vs lineages (columns)
tissue_lineage_counts <- table(tm.facs$tissue, tm.facs$cell_lineage)
print(tissue_lineage_counts)

# ________________
# Quality control
# ________________

# Calculate QC metrics
tm.facs <- addPerCellQC(tm.facs, 
                        subsets=list(mito=grep("^mt-", rownames(tm.facs))))

# Identify outliers based on library size, gene count and mitochondrial percentage
qc_lib <- isOutlier(tm.facs$sum, log=TRUE, type="lower", nmads = 3)
qc_nexprs <- isOutlier(tm.facs$detected, log=TRUE, type="lower", nmads = 3)
qc_mito <- isOutlier(tm.facs$subsets_mito_percent, type="higher", nmads = 3)

# Flag outlier cells
tm.facs$discard <- qc_lib | qc_nexprs | qc_mito

# Summarize number of cells removed by QC
message("Number of cells flagged as outliers: ", sum(tm.facs$discard), 
        " out of ", ncol(tm.facs), " (", 
        round(sum(tm.facs$discard)/ncol(tm.facs)*100, 2), "%)")

# Remove low-quality cells
tm.facs <- tm.facs[, !tm.facs$discard]
message("Dimensions after QC: ", paste(dim(tm.facs), collapse=" x "))

# ______________
# Normalization
# ______________

# Log-normalization
tm.facs <- logNormCounts(tm.facs)
message("Normalization complete. assays available: ", 
        paste(assayNames(tm.facs), collapse=", "))

# __________________________________
# Extract myeloid cells from Marrow
# __________________________________

myeloid_cell_types <- c("macrophage", "alveolar macrophage", "monocyte", 
                        "non-classical monocyte", "classical monocyte", 
                        "dendritic cell", "myeloid cell", "granulocyte",
                        "Langerhans cell", "leukocyte", "basophil", 
                        "mast cell", "blood cell", 
                        "granulocytopoietic cell", "promonocyte")

# Focus on Marrow tissue only
marrow_cells <- tm.facs[, tm.facs$tissue == "Marrow"]
rm(tm.facs)
marrow_myeloid <- marrow_cells[, marrow_cells$cell_ontology_class %in% myeloid_cell_types]

# Check number of cells per type
celltype_counts <- table(marrow_myeloid$cell_ontology_class)
print(celltype_counts)

# ________________________________________
# Split into reference and query datasets
# ________________________________________

cell_indices <- 1:ncol(marrow_myeloid)
reference_indices <- sample(cell_indices, size = floor(length(cell_indices) * 0.7))
query_indices <- setdiff(cell_indices, reference_indices)

reference_cells <- marrow_myeloid[, reference_indices]
query_cells <- marrow_myeloid[, query_indices]
reference_cells_subset <- reference_cells[, reference_cells$cell_ontology_class != "promonocyte"]

message("Reference dataset: ", ncol(reference_cells), " cells")
message("Reference subset dataset: ", ncol(reference_cells_subset), " cells")
message("Query dataset: ", ncol(query_cells), " cells")

# Check distribution in reference and query
ref_dist <- table(reference_cells$cell_ontology_class)
ref_subset_dist <- table(reference_cells_subset$cell_ontology_class)
query_dist <- table(query_cells$cell_ontology_class)
print("Reference distribution:")
print(ref_dist)
print("Reference subset distribution:")
print(ref_subset_dist)
print("Query distribution:")
print(query_dist)

# __________________
# Feature selection
# __________________

# Identify variable genes - full data
var_decomp <- modelGeneVar(marrow_myeloid)
top_hvgs <- getTopHVGs(var_decomp, n = 2000)
message("Number of HVGs selected: ", length(top_hvgs))

# Identify variable genes - subset data
var_decomp_subset <- modelGeneVar(marrow_myeloid[, marrow_myeloid$cell_ontology_class != "promonocyte"])
top_hvgs_subset <- getTopHVGs(var_decomp_subset, n = 2000)
message("Number of HVGs selected for subset: ", length(top_hvgs_subset))

# ____________________
# Dimension reduction
# ____________________

# PCA
reference_cells <- runPCA(reference_cells, subset_row = top_hvgs)
reference_cells_subset <- runPCA(reference_cells_subset, subset_row = top_hvgs_subset)

# UMAP
reference_cells <- runUMAP(reference_cells, dimred = "PCA", n_dimred = 30)
umap_plot <- plotReducedDim(reference_cells, 
                            dimred = "UMAP", 
                            colour_by = "cell_ontology_class")
print(umap_plot)
reference_cells_subset <- runUMAP(reference_cells_subset, dimred = "PCA", n_dimred = 30)
umap_plot_subset <- plotReducedDim(reference_cells_subset, 
                                   dimred = "UMAP", 
                                   colour_by = "cell_ontology_class")
print(umap_plot_subset)

# ___________________________
# Perform SingleR annotation
# ___________________________

# Use HVGs for SingleR
common_hvgs <- intersect(top_hvgs, rownames(marrow_myeloid))
common_hvgs_subset <- intersect(top_hvgs_subset, rownames(marrow_myeloid))

# Run SingleR - full data
start_time <- Sys.time()
results <- SingleR(
    test = query_cells[common_hvgs,],
    ref = reference_cells[common_hvgs,],
    labels = reference_cells$cell_ontology_class
)
end_time <- Sys.time()
message("SingleR annotation completed in ", round(difftime(end_time, start_time, units="mins"), 2), " minutes")

# Add SingleR labels to query cells
query_cells$SingleR_labels <- results$labels

# Calculate scores per cell type
score_heatmap <- plotScoreHeatmap(results)
print(score_heatmap)

# Run SingleR - subset of data
start_time <- Sys.time()
results_subset <- SingleR(
    test = query_cells[common_hvgs_subset,],
    ref = reference_cells[common_hvgs_subset, 
                          reference_cells$cell_ontology_class != "promonocyte"],
    labels = reference_cells$cell_ontology_class[
        reference_cells$cell_ontology_class != "promonocyte"]
)
end_time <- Sys.time()
message("SingleR annotation completed in ", round(difftime(end_time, start_time, units="mins"), 2), " minutes")

# Add SingleR labels to query cells
query_cells$SingleR_labels_subset <- results_subset$labels

# Calculate scores per cell type
score_heatmap_subset <- plotScoreHeatmap(results_subset)
print(score_heatmap_subset)

# _____________________________
# Evaluate annotation accuracy
# _____________________________

# Compare original labels with predicted labels
confusion_matrix <- table(Original = query_cells$cell_ontology_class, 
                          Predicted = query_cells$SingleR_labels_subset)
print(confusion_matrix)

# Calculate overall accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
message("Overall annotation accuracy: ", round(accuracy * 100, 2), "%")

# Visualize agreement/disagreement
query_cells$label_match <- query_cells$cell_ontology_class == query_cells$SingleR_labels_subset

# ______________________________________________
# Save preprocessed data for scDiagnostics demo
# ______________________________________________

# Save the processed objects for your package demonstration
saveRDS(reference_cells, "data/reference_marrow_myeloid.rds")
saveRDS(reference_cells_subset, "data/reference_subset_marrow_myeloid.rds")
saveRDS(query_cells, "data/query_marrow_myeloid.rds")

message("Data processing complete. Objects saved for scdiagnostics demonstration.")







