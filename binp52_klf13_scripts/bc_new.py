import os

from sklearn.decomposition import TruncatedSVD
from sklearn.kernel_approximation import svd
os.environ['RPY2_CFFI_MODE'] = 'ABI'  # harmless even if you don't use R

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import anndata
from sklearn.metrics import silhouette_score
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------
# Settings
# --------------------------
INPUT_FILE = "combined_8_samples_all_layers_outer.h5ad"
NORM_LAYER = "analytic_pearson_residuals"   # single layer to correct
LEIDEN_RESOLUTIONS = [0.25, 0.5, 0.7, 1, 2]

# --------------------------
# Helpers
# --------------------------
def load_for_combat(input_file: str, norm_layer_name: str) -> anndata.AnnData:
    print(f"Loading 10K HVGs from '{norm_layer_name}' (ultra low mem)...")
    
    # Load metadata only
    adata_full = sc.read_h5ad(input_file, backed='r')
    
    # **STRATEGY 1**: Top 10K genes by mean (perfect proxy for HVGs, no computation needed)
    print("Selecting top 10K genes by expression (fast HVG proxy)...")
    layer_means = np.array(adata_full.layers[norm_layer_name].mean(axis=0)).flatten()
    hvg_idx = np.argsort(layer_means)[-10000:]
    hvg_genes = adata_full.var_names[hvg_idx]
    
    # **STRATEGY 2**: Load ONLY these 10K genes directly
    print("Loading 10K genes only...")
    X_hvg = adata_full.layers[norm_layer_name][:, hvg_idx]
    
    adata = anndata.AnnData(
        X=X_hvg,
        obs=adata_full.obs.copy(),
        var=adata_full.var.loc[hvg_genes].copy(),
        obsm=adata_full.obsm.copy() if 'obsm' in adata_full.__dict__ else {}
    )
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes (~150MB)")
    
    adata_full.file.close()
    return adata


def run_leiden_silhouette(adata: anndata.AnnData,
                          embedding_key: str,
                          label_prefix: str):
    """Run Leiden at multiple resolutions and compute silhouette scores."""
    sil_scores = {}
    for res in LEIDEN_RESOLUTIONS:
        r_str = str(res).replace('.', '_')
        leiden_key = f"{label_prefix}_{r_str}"
        sc.tl.leiden(adata, resolution=res, key_added=leiden_key)
        labels = adata.obs[leiden_key].astype(str)
        score = silhouette_score(adata.obsm[embedding_key], labels)
        sil_scores[res] = score
        print(f"Leiden res={res}: silhouette={score:.4f}")

    best_res = max(sil_scores, key=sil_scores.get)
    best_key = f"{label_prefix}_{str(best_res).replace('.', '_')}"
    print(f"Best resolution ({label_prefix}): {best_res}")
    return best_res, best_key, sil_scores


def run_combat_pipeline(adata: anndata.AnnData,
                        norm_name: str):
    """
    ComBat batch correction + metrics + before/after UMAP on 3000 HVGs.
    Uses scanpy's ComBat implementation.
    """
    print("\n" + "=" * 60)
    print(f"=== ComBat batch correction on {norm_name} ===")
    print("=" * 60)

    # PCA / neighbors / UMAP BEFORE ComBat
    print("Running PCA on HVGs...")
    #sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
    
    from sklearn.decomposition import TruncatedSVD
    svd = TruncatedSVD(n_components=50, algorithm='arpack')
    adata.obsm['X_pca'] = svd.fit_transform(adata.X)
    print(f"TruncatedSVD complete: {svd.explained_variance_ratio_.sum():.3f} variance explained")

    print("Computing neighbors/UMAP BEFORE ComBat...")
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=0)
    adata.obsm["X_umap_before"] = adata.obsm["X_umap"].copy()

    # ComBat expects dense arrays in Scanpy implementation[web:46]
    print("Running ComBat (dense X)...")
    if sp.issparse(adata.X):
        adata.X = adata.X.toarray()
    sc.pp.combat(adata, key="batch")

    # PCA / neighbors / UMAP AFTER ComBat
    print("Computing PCA/neighbors/UMAP AFTER ComBat...")
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=0)
    adata.obsm["X_umap_after"] = adata.obsm["X_umap"].copy()

    # Leiden + silhouette on corrected PCA embedding
    print("Running Leiden + silhouette on ComBat PCs...")
    best_res, best_leiden_key, sil_scores = run_leiden_silhouette(
        adata,
        embedding_key="X_pca",
        label_prefix=f"leiden_combat_{norm_name}",
    )

    # Save corrected object
    out_file = f"combined_8_samples_batch_corrected_combat_{norm_name}.h5ad"
    adata.write(out_file)
    print(f"Saved ComBat-corrected data: {out_file}")

    print("=" * 60 + "\n")

    return {
        "norm_layer": norm_name,
        "best_res": best_res,
        "best_leiden_key": best_leiden_key,
        "sil_scores": sil_scores,
        "outfile": out_file,
    }


# --------------------------
# Main (single norm layer)
# --------------------------
if __name__ == "__main__":
    print("Available layers in input file:")
    adata_test = sc.read_h5ad(INPUT_FILE, backed='r')
    print(list(adata_test.layers.keys()))
    adata_test.file.close()

    # Load + prepare for ONE norm layer with HVGs
    print(f"Processing normalization layer: {NORM_LAYER}")
    adata = load_for_combat(INPUT_FILE, NORM_LAYER)

    # Run ComBat pipeline
    stats = run_combat_pipeline(adata, NORM_LAYER)

    # BEFORE ComBat
    svd_before = TruncatedSVD(n_components=50, algorithm='arpack')
    adata.obsm['X_pca'] = svd_before.fit_transform(adata.X)

    # AFTER ComBat  
    svd_after = TruncatedSVD(n_components=50, algorithm='arpack') 
    adata.obsm['X_pca'] = svd_after.fit_transform(adata.X)

    print("Done.")
