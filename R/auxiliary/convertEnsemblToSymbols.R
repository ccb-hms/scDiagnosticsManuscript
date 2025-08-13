# ---------------------------------------------------
# Convert Ensembl IDs to Gene Symbols for SCE Object
# ---------------------------------------------------

# Required libraries
library(biomaRt)             
library(SingleCellExperiment)

# ---------------------------------------------------

# Convert Ensembl IDs to gene symbols
convertEnsemblToSymbols <- function(sce_obj) {
    cat("Converting Ensembl IDs to gene symbols...\n")
    
    ensembl_ids <- rownames(sce_obj)
    
    # Check if IDs look like Ensembl IDs
    if(!any(grepl("^ENS", ensembl_ids))) {
        cat("Warning: Row names don't appear to be Ensembl IDs\n")
        return(sce_obj)
    }
    
    # Connect to biomaRt with error handling
    mart <- tryCatch({
        useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }, error = function(e) {
        cat("Error connecting to biomaRt:", e$message, "\n")
        return(NULL)
    })
    
    if(is.null(mart)) {
        cat("Failed to connect to biomaRt, returning original object\n")
        return(sce_obj)
    }
        
    # Get gene symbols
    gene_info <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id", 
        values = ensembl_ids,
        mart = mart
    )
    
    # Match to our genes
    matched_symbols <- gene_info$hgnc_symbol[match(ensembl_ids, gene_info$ensembl_gene_id)]
    
    # Use Ensembl ID if no symbol found
    gene_symbols <- ifelse(is.na(matched_symbols) | matched_symbols == "", 
                           ensembl_ids, matched_symbols)
    
    # Handle duplicated symbols by making them unique
    gene_symbols <- make.unique(gene_symbols, sep = "_")
    
    cat("Converted", sum(!is.na(matched_symbols) & matched_symbols != ""), 
        "Ensembl IDs to gene symbols\n")
    
    # Update rownames and rowData
    rownames(sce_obj) <- gene_symbols
    rowData(sce_obj)$ensembl_id <- ensembl_ids
    rowData(sce_obj)$gene_symbol <- gene_symbols
    
    return(sce_obj)
}