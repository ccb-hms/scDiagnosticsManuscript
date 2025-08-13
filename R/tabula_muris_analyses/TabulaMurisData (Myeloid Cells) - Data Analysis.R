# -----------------------------------------------
# Process TabulaMurisData - Marrow Myeloid Cells
# -----------------------------------------------

# Load libraries
library(scDiagnostics)

# Source files


# -----------------------------------------------

# Load reference and query SCE objects
query_cells <- readRDS("data/query_marrow_myeloid.rds")
reference_cells <- readRDS("data/reference_marrow_myeloid.rds")
reference_cells_subset <- readRDS("data/reference_subset_marrow_myeloid.rds")

# PCA diagnostics
plotCellTypePCA(query_data = query_cells,
                reference_data = reference_cells, 
                query_cell_type_col = "SingleR_labels_subset",
                ref_cell_type_col = "cell_ontology_class", 
                pc_subset = 1:3,
                assay_name = "logcounts", 
                diagonal_facet = "ridge", 
                upper_facet = "blank")

# Discriminant diagnostics
disc_output <- calculateDiscriminantSpace(query_data = query_cells,
                                          reference_data = reference_cells_subset, 
                                          query_cell_type_col = "SingleR_labels_subset",
                                          ref_cell_type_col = "cell_ontology_class", 
                                          n_tree = 500,
                                          n_top = 50,
                                          eigen_threshold  = 1e-1,
                                          calculate_metrics = FALSE,
                                          alpha = 0.01)
plot(disc_output,
     diagonal_facet = "ridge", 
     upper_facet = "blank", 
     cell_types = NULL)

# SIR diagnostics
sir_output <- calculateSIRSpace(query_data = query_cells,
                                reference_data = reference_cells_subset, 
                                query_cell_type_col = "SingleR_labels_subset",
                                ref_cell_type_col = "cell_ontology_class", 
                                multiple_cond_means = TRUE,
                                cumulative_variance_threshold = 0.9,
                                n_neighbor = 1)
plot(sir_output,
     sir_subset = 6:10,
     cell_types = NULL,
     lower_facet = "scatter",
     diagonal_facet = "ridge",
     upper_facet = "blank")

# Pairwise distances diagnostics
plotPairwiseDistancesDensity(query_data = query_cells,
                             reference_data = reference_cells_subset, 
                             query_cell_type_col = "SingleR_labels_subset",
                             ref_cell_type_col = "cell_ontology_class", 
                             pc_subset = 1:5, 
                             cell_type_query = "monocyte",
                             cell_type_ref = "monocyte",
                             distance_metric = "euclidean", 
                             bandwidth = 0.45)

# Pairwise correlation diagnostics
cor_plot <- 
    calculateAveragePairwiseCorrelation(query_data = query_cells,
                                        reference_data = reference_cells_subset, 
                                        query_cell_type_col = "SingleR_labels_subset",
                                        ref_cell_type_col = "cell_ontology_class", 
                                        pc_subset = 1:10, 
                                        correlation_method = "pearson")
plot(cor_plot)

# Isolation forest diagnostics
anomaly_output <- detectAnomaly(query_data = query_cells,
                                reference_data = reference_cells_subset, 
                                query_cell_type_col = "SingleR_labels_subset",
                                ref_cell_type_col = "cell_ontology_class", 
                                pc_subset = 1:10,
                                n_tree = 500,
                                anomaly_threshold = 0.525)
plot(anomaly_output,
     cell_type = "monocyte",
     pc_subset = 1:5,
     data_type = "query",
     diagonal_facet = "density", 
     upper_facet = "blank")

# Cell distances disgnostics
distance_data <- calculateCellDistances(query_data = query_cells,
                                        reference_data = reference_cells_subset, 
                                        query_cell_type_col = "SingleR_labels_subset",
                                        ref_cell_type_col = "cell_ontology_class", 
                                        pc_subset = 1:10)
monocyte_top6_anomalies <- 
    names(sort(anomaly_output$monocyte$query_anomaly_scores, decreasing = TRUE)[1:6])
plot(distance_data, ref_cell_type = "monocyte", cell_names = monocyte_top6_anomalies)

