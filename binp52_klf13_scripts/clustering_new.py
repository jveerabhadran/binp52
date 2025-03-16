import scanpy as sc
from sklearn.metrics import silhouette_score
import numpy as np

# Set Scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)


# Load the preprocessed AnnData object
adata = sc.read("dimensionality_reduced_scran_combat.h5ad")


# Step 10.1: Construct the KNN graph on the PCA-reduced space
print("Constructing KNN graph...")
sc.pp.neighbors(adata, n_pcs=30) 


# Step 10.2: Compute UMAP embedding for visualization
print("Computing UMAP embedding...")
sc.tl.umap(adata)


# Step 10.3: Run Leiden clustering at multiple resolutions
print("Running Leiden clustering at multiple resolutions...")
sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
sc.tl.leiden(adata, key_added="leiden_res2", resolution=2.0)


# Step 10.4: Visualize clustering results with UMAP
print("Visualizing clustering results...")
sc.pl.umap(
    adata,
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1", "leiden_res2"],
    legend_loc="on data",
)

res_keys = ["leiden_res0_25", "leiden_res0_5", "leiden_res1", "leiden_res2"]

for key in res_keys:
    # Silhouette needs clusters as integers, so convert Leiden labels if needed
    labels = adata.obs[key].astype(int)
    X = adata.obsm["X_pca"]
    if len(np.unique(labels)) > 1 and len(np.unique(labels)) < len(labels):  # Need >1 cluster & <N clusters
        sil_score = silhouette_score(X, labels)
        print(f"Silhouette score for {key}: {sil_score:.3f}")
    else:
        print(f"{key} not valid for silhouette analysis (all same or all unique).")


# Save the updated AnnData object with clustering results
output_file = "ko_scran_clustered.h5ad"
adata.write(output_file)


print(f"Processing completed. Saved output to {output_file}.")
