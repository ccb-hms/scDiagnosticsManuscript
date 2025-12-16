# -----------------------------------------------------------
# MERFISH Mouse Colon IBD - CellTypist Cell Type Annotation
# -----------------------------------------------------------

# Load libraries
import pandas as pd
import scanpy as sc
import celltypist
import numpy as np
from datetime import datetime
import os

# -----------------------------------------------------------

# Configuration
# Ensuring paths match the output from the R Data Prep script
ANNDATA_REF_FILE = "data/merfish/reference_data.h5ad"
ANNDATA_QUERY_FILE = "data/merfish/query_data.h5ad" 
ANNOTATED_QUERY_FILE = "data/merfish/annotated_query_data.h5ad"
CELL_TYPE_KEY = "tier2"  

# Set scanpy settings
sc.settings.verbosity = 2  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_and_preprocess_data(adata, target_sum=1e4):
    """Load and preprocess AnnData object"""
    print("Data shape: {adata.shape}")
    
    # Check if data is sparse or dense
    if hasattr(adata.X, "toarray"):
        print("Data type: Sparse matrix")
    else:
        print("Data type: Dense matrix")
    
    # Basic preprocessing
    # Note: MERFISH data is often already normalized, but for CellTypist 
    # ensuring standard log1p normalization is recommended if starting from counts.
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    return adata

def train_celltypist_model(adata_ref, cell_type_key=CELL_TYPE_KEY):
    """Train CellTypist model from reference data"""
    print("Training CellTypist model...")
    print(f"Cell types in reference: {adata_ref.obs[cell_type_key].value_counts()}")
    
    start_time = datetime.now()
    
    # Train the model
    # Note: Since MERFISH panel has fewer genes than whole transcriptome, 
    # we train a custom model specifically on the available genes.
    model = celltypist.train(
        adata_ref,
        labels=cell_type_key,
        n_jobs=-1,  # Use all available cores
        max_iter=100,
        mini_batch=True  # For efficiency
    )
    
    end_time = datetime.now()
    print(f"Training completed in: {end_time - start_time}")
    
    return model

def annotate_query_data(model, adata_query):
    """Annotate query data using trained model"""
    print("Annotating query data...")
    start_time = datetime.now()
    
    # Predict cell types 
    predictions = celltypist.annotate(
        adata_query,
        model=model,
        majority_voting=False  # Disable to avoid memory issues and ensure direct mapping
    )
    
    end_time = datetime.now()
    print(f"Annotation completed in: {end_time - start_time}")
    
    # Get annotated data
    adata_annotated = predictions.to_adata()
    
    # Print prediction summary
    print("\nPrediction summary:")
    if 'predicted_labels' in adata_annotated.obs:
        print(f"Predicted labels: {adata_annotated.obs['predicted_labels'].value_counts().head()}")
    
    # Print confidence score statistics
    if 'conf_score' in adata_annotated.obs.columns:
        conf_scores = adata_annotated.obs['conf_score']
        print(f"\nConfidence scores - Mean: {conf_scores.mean():.3f}, Std: {conf_scores.std():.3f}")
        print(f"Low confidence cells (<0.5): {sum(conf_scores < 0.5)} ({100*sum(conf_scores < 0.5)/len(conf_scores):.1f}%)")
    
    return adata_annotated
        
def main():
    """Main annotation workflow"""
    print("=== MERFISH CellTypist Annotation Workflow ===")
    total_start = datetime.now()
    
    # Verify files exist
    if not os.path.exists(ANNDATA_REF_FILE):
        raise FileNotFoundError(f"Reference file not found: {ANNDATA_REF_FILE}")
    if not os.path.exists(ANNDATA_QUERY_FILE):
        raise FileNotFoundError(f"Query file not found: {ANNDATA_QUERY_FILE}")
    
    # Load data
    print("Loading reference data...")
    adata_ref = sc.read_h5ad(ANNDATA_REF_FILE)
    
    print("Loading query data...")
    adata_query = sc.read_h5ad(ANNDATA_QUERY_FILE)
    
    # Preprocess data
    print("\nPreprocessing reference data...")
    adata_ref = load_and_preprocess_data(adata_ref)
    
    print("Preprocessing query data...")
    adata_query = load_and_preprocess_data(adata_query)
    
    # Train model
    print("\n" + "="*50)
    model = train_celltypist_model(adata_ref)
    
    # Annotate query
    print("\n" + "="*50) 
    adata_annotated = annotate_query_data(model, adata_query)
    
    # Save results
    print(f"\nSaving annotated data to {ANNOTATED_QUERY_FILE}...")
    adata_annotated.write_h5ad(ANNOTATED_QUERY_FILE)
    
    total_end = datetime.now()
    print(f"\n=== Total workflow time: {total_end - total_start} ===")
    print("Annotation complete!")

if __name__ == "__main__":
    main()