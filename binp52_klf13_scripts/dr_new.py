import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# Set Scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# Load the preprocessed AnnData object
adata = sc.read(
    filename="Combined\combined_feature_selected.h5ad",
)

# Use the normalized representation for dimensionality reduction
adata.X = adata.layers["log1p_norm"]

# Step 1: PCA
print("Running PCA...")
adata.var["highly_variable"] = adata.var["highly_deviant"]  # Use highly deviant genes for PCA
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)  # Perform PCA
sc.pl.pca_scatter(adata, color="total_counts")  # Visualize PCA scatter plot

# Step 2: t-SNE
print("Running t-SNE...")
sc.tl.tsne(adata, use_rep="X_pca")  # Perform t-SNE using PCA representation
sc.pl.tsne(adata, color="total_counts")  # Visualize t-SNE scatter plot

# Step 3: UMAP
print("Running UMAP...")
sc.pp.neighbors(adata)  # Compute neighborhood graph based on PCA
sc.tl.umap(adata)  # Perform UMAP
sc.pl.umap(adata, color="total_counts")  # Visualize UMAP scatter plot

# Step 4: Inspecting Quality Control Metrics in UMAP
print("Inspecting quality control metrics...")
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
)  # Visualize QC metrics on UMAP

# Save the updated AnnData object with dimensionality reduction results
output_file = "Combined\combined_dimensionality_reduced.h5ad"
adata.write(output_file)

print(f"Processing completed. Saved output to {output_file}.")