# ---------------------------------------------------------
# MERFISH Mouse Colon - scVI/scArches Annotation and UMAP
# ---------------------------------------------------------

import scvi
import scarches as sca
import scanpy as sc
import anndata
import pandas as pd

# --- Configuration ---
ANNDATA_REF_FILE = "data/merfish/scvi_reference_data.h5ad"
ANNDATA_QUERY_FILE = "data/merfish/scvi_query_data.h5ad"
QUERY_PRED_CSV_FILE = "data/merfish/scvi_predictions.csv"
REF_UMAP_CSV_FILE = "data/merfish/scvi_reference_umap.csv"
CELL_TYPE_KEY = "tier2"
BATCH_KEY = "sample_id"

def main():
    print("--- Starting scVI/scArches Annotation Workflow (MERFISH) ---")
    scvi.settings.seed = 0

    # 1. Load Data
    adata_ref = sc.read_h5ad(ANNDATA_REF_FILE)
    adata_query = sc.read_h5ad(ANNDATA_QUERY_FILE)
    
    original_ref_barcodes = adata_ref.obs.index.tolist()
    original_query_barcodes = adata_query.obs.index.tolist()
    reference_labels = adata_ref.obs[CELL_TYPE_KEY].copy()

    # 2. Setup and Train Reference Model
    print("\n--- Training Reference scVI Model ---")
    scvi.model.SCVI.setup_anndata(adata_ref, batch_key=BATCH_KEY)
    vae_ref = scvi.model.SCVI(adata_ref)
    vae_ref.train()

    # 3. Map Query Data
    print("\n--- Mapping Query Data with scArches ---")
    vae_query = scvi.model.SCVI.load_query_data(adata_query, vae_ref)
    vae_query.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))

    # 4. Get Latent Representations and Compute Joint UMAP
    print("\n--- Generating Latent Space and Joint UMAP ---")
    adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
    adata_query.obsm["X_scVI"] = vae_query.get_latent_representation()

    adata_full = anndata.concat(
        {"ref": adata_ref, "query": adata_query},
        label="data_source", index_unique=None
    )
    
    sc.pp.neighbors(adata_full, use_rep="X_scVI")
    sc.tl.umap(adata_full)

    # 5. Separate and perform annotation with scArches
    print("\n--- Annotating Query Cells via scArches weighted k-NN ---")
    adata_ref_processed = adata_full[adata_full.obs.data_source == 'ref'].copy()
    adata_query_processed = adata_full[adata_full.obs.data_source == 'query'].copy()

    # Add cell type labels to reference
    adata_ref_processed.obs[CELL_TYPE_KEY] = reference_labels

    # Use scArches weighted k-NN trainer
    knn_transformer = sca.utils.knn.weighted_knn_trainer(
        train_adata=adata_ref_processed,
        train_adata_emb="X_scVI",
        n_neighbors=50,
    )

    # Transfer labels with uncertainty
    labels, uncert = sca.utils.knn.weighted_knn_transfer(
        query_adata=adata_query_processed,
        query_adata_emb="X_scVI",
        label_keys=CELL_TYPE_KEY,
        knn_model=knn_transformer,
        ref_adata_obs=adata_ref_processed.obs,
    )

    query_pred = labels.loc[adata_query_processed.obs.index, CELL_TYPE_KEY].values
    query_confidence = uncert.loc[adata_query_processed.obs.index, CELL_TYPE_KEY].values

    # 6. Create pandas DataFrames
    print("--- Assembling results into DataFrames ---")
    
    query_results_df = pd.DataFrame({
        'barcode': adata_query_processed.obs.index,
        'scvi_prediction': query_pred,
        'scvi_confidence': query_confidence,
        'UMAP_scVI_1': adata_query_processed.obsm["X_umap"][:, 0],
        'UMAP_scVI_2': adata_query_processed.obsm["X_umap"][:, 1]
    }).set_index('barcode').loc[original_query_barcodes]

    ref_results_df = pd.DataFrame({
        'barcode': adata_ref_processed.obs.index,
        'UMAP_scVI_1': adata_ref_processed.obsm["X_umap"][:, 0],
        'UMAP_scVI_2': adata_ref_processed.obsm["X_umap"][:, 1]
    }).set_index('barcode').loc[original_ref_barcodes]

    # 7. Save results
    print(f"\nSaving query predictions and UMAP to {QUERY_PRED_CSV_FILE}...")
    query_results_df.to_csv(QUERY_PRED_CSV_FILE)
    
    print(f"Saving reference UMAP to {REF_UMAP_CSV_FILE}...")
    ref_results_df.to_csv(REF_UMAP_CSV_FILE)
    
    print("✓ Python workflow complete.")

if __name__ == "__main__":
    main()