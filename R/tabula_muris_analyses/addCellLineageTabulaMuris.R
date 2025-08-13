# -----------------------------------------------------
# Add Cell Lineage Classification - Tabula Muris data
# -----------------------------------------------------


# This function adds a 'cell_lineage' column to tm.facs that classifies each cell
# into one of four broader categories
addCellLineageTabulaMuris <- function(sce) {
    
    # Create a new column with default value of "other"
    sce$cell_lineage <- "other"
    
    # Define the cell types for each lineage
    myeloid_cell_types <- c(
        "macrophage", "alveolar macrophage", "monocyte", 
        "non-classical monocyte", "classical monocyte", 
        "dendritic cell", "myeloid cell", "granulocyte",
        "Langerhans cell", "leukocyte", "basophil", 
        "mast cell", "blood cell", 
        "granulocytopoietic cell", "promonocyte" 
    )
    
    # New category for hematopoietic progenitors
    hematopoietic_progenitor_types <- c(
        "hematopoietic precursor cell", 
        "erythroblast", 
        "proerythroblast"
    )
    
    b_cell_types <- c(
        "Fraction A pre-pro B cell", "early pro-B cell", 
        "late pro-B cell", "immature B cell", "B cell"
    )
    
    t_cell_types <- c(
        "immature T cell", "DN1 thymic pro-T cell", 
        "T cell", "natural killer cell"
    )
    
    epithelial_cell_types <- c(
        "bladder urothelial cell", "kidney collecting duct epithelial cell",
        "kidney proximal straight tubule epithelial cell", 
        "kidney loop of Henle ascending limb epithelial cell",
        "duct epithelial cell", "epithelial cell", "basal cell of epidermis",
        "keratinocyte", "luminal epithelial cell of mammary gland",
        "type II pneumocyte", "ciliated columnar cell of tracheobronchial tree"
    )
    
    # Assign lineage based on cell type
    sce$cell_lineage[sce$cell_ontology_class %in% myeloid_cell_types] <- "myeloid_cell_lineage"
    sce$cell_lineage[sce$cell_ontology_class %in% hematopoietic_progenitor_types] <- "hematopoietic_progenitors"
    sce$cell_lineage[sce$cell_ontology_class %in% b_cell_types] <- "b_cell_development_trajectory"
    sce$cell_lineage[sce$cell_ontology_class %in% t_cell_types] <- "t_cell_subtype"
    sce$cell_lineage[sce$cell_ontology_class %in% epithelial_cell_types] <- "epithelial_cell_type"
    
    # Return the modified SingleCellExperiment
    return(sce)
}
