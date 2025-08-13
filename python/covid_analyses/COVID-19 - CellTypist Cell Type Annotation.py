# ----------------------------------------------
# COVID-19 - CellTypist Cell Type Annotation
# ----------------------------------------------

# Load libraries
import pandas as pd
import scanpy as sc
import celltypist
import numpy as np
from datetime import datetime

# ----------------------------------------------

# Configuration
ANNDATA_REF_FILE = "data/reference_data.h5ad"
ANNDATA_QUERY_FILE = "data/query_data.h5ad" 
ANNOTATED_QUERY_FILE = "data/annotated_query_data.h5ad"
CELL_TYPE_KEY = "cell_type"

# Set scanpy settings
sc.settings.verbosity = 2  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_and_preprocess_data(adata, target_sum=1e4):
    """Load and preprocess AnnData object"""
    print(f"Data shape: {adata.shape}")
    print(f"Data type: {type(adata.X)}")
    
    # Basic preprocessing
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    return adata

def train_celltypist_model(adata_ref, cell_type_key=CELL_TYPE_KEY):
    """Train CellTypist model from reference data"""
    print("Training CellTypist model...")
    print(f"Cell types in reference: {adata_ref.obs[cell_type_key].value_counts()}")
    
    start_time = datetime.now()
    
    # Train the model
    model = celltypist.train(
        adata_ref,
        labels=cell_type_key,
        n_jobs=-1,  # Use all available cores
        max_iter=100,
        mini_batch=True  # For efficiency with large datasets
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
        majority_voting=True,  # More robust predictions
        over_clustering='auto'  # Automatic over-clustering
    )
    
    end_time = datetime.now()
    print(f"Annotation completed in: {end_time - start_time}")
    
    # Get annotated data
    adata_annotated = predictions.to_adata()
    
    # Print prediction summary
    print("\nPrediction summary:")
    print(f"Predicted labels: {adata_annotated.obs['predicted_labels'].value_counts()}")
    
    if 'majority_voting' in adata_annotated.obs.columns:
        print(f"Majority voting labels: {adata_annotated.obs['majority_voting'].value_counts()}")
    
    # Print confidence score statistics
    if 'conf_score' in adata_annotated.obs.columns:
        conf_scores = adata_annotated.obs['conf_score']
        print(f"\nConfidence scores - Mean: {conf_scores.mean():.3f}, Std: {conf_scores.std():.3f}")
        print(f"Low confidence cells (<0.5): {sum(conf_scores < 0.5)} ({100*sum(conf_scores < 0.5)/len(conf_scores):.1f}%)")
    
    return adata_annotated

def main():
    """Main annotation workflow"""
    print("=== COVID-19 CellTypist Annotation Workflow ===")
    total_start = datetime.now()
    
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