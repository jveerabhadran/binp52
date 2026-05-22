import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

# UPDATE: Your actual file
adata = sc.read_h5ad("combined_8_samples_batch_corrected_combat_log1p_norm.h5ad")

# Use scran layer (matches your correction)
#adata.X = adata.layers["scran_normalization"].copy()  # Fixed .copy()

# FIXED - ComBat already corrected .X, no layer copy needed
print("Using ComBat-corrected data in adata.X (no layer copy needed)")
print(f"Shape: {adata.X.shape}")

# PCA using existing HVGs
if "highly_deviant" not in adata.var:
    adata.var["highly_deviant"] = adata.var_names.isin(adata.var_names)  # All genes
sc.pp.pca(adata, svd_solver="arpack")
# FIXED plot
#sc.pl.pca_variance_explained(adata, log=True, n_pcs=50)  # Elbow plot

# Original PCA scatter (optional)
sc.pl.pca_scatter(adata, color="total_counts", components=['1,2'])


# t-SNE
sc.tl.tsne(adata, use_rep="X_pca")
sc.pl.tsne(adata, color="total_counts")

# UMAP (already exists from ComBat, but recompute)
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)
sc.pl.umap(adata, color="total_counts")

# QC metrics
sc.pl.umap(adata, color=["total_counts", "pct_counts_mt", "scDblFinder_class"], ncols=2)

# Save
output_file = "lognorm_dimensionality_reduced.h5ad"
adata.write(output_file)
print(f"Saved: {output_file}")

