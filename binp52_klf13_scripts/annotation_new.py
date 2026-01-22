import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# Set figure parameters
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)


# Load the AnnData object
adata = sc.read_h5ad("lognorm_clustered.h5ad")

print("Layers:", list(adata.layers.keys()))
print("Shape:", adata.n_obs, "cells x", adata.n_vars, "genes")

# Store raw counts
adata.layers["counts"] = adata.X.copy()  # Ensure to copy to avoid modifying the original data

# Set APR residuals for downstream annotation/plotting
#adata.X = adata.layers["analytic_pearson_residuals"].copy()
#layer_for_plot = "analytic_pearson_residuals"

adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=10.0, neginf=0.0)
layer_for_plot = None  # Fix: No layers → use adata.X


# Step 1: Clustering Analysis
# Perform Leiden clustering at different resolutions and visualize results in UMAP
#sc.tl.leiden(adata, resolution=1, key_added="leiden_1")
#sc.pl.umap(adata, color="leiden_1", title="Leiden Clustering (Resolution 1)")

print("\nRunning Leiden clustering at resolution=0.25...")
sc.pp.neighbors(adata, use_rep='X_umap')  # Use existing UMAP
sc.tl.leiden(adata, resolution=0.25, key_added="leiden_res0_25")
print("Leiden res0_25 clusters:", adata.obs['leiden_res0_25'].nunique())

# Plot the new clustering
sc.pl.umap(adata, color="leiden_res0_25", legend_loc='on data', frameon=False, title="Leiden Clustering (Res 0.25)")
# Step 2: Cell Cycle Marker Analysis
# Define marker genes for S phase and G2/M phase
s_genes = [
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6",
    "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1",
    "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2"
]

g2m_genes = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2",
    "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Fam64a", "Smc4", "CcnB2", "Ckap2l", "Ckap2"
]

# Use lowercase (mouse genes)
s_genes_mouse = [gene.lower() for gene in s_genes]
g2m_genes_mouse = [gene.lower() for gene in g2m_genes]

s_genes_filtered = [gene for gene in s_genes_mouse if gene in adata.var.index]
g2m_genes_filtered = [gene for gene in g2m_genes_mouse if gene in adata.var.index]

print("Filtered S genes:", s_genes_filtered)
print("Filtered G2M genes:", g2m_genes_filtered)

# Skip cell cycle completely - your 3k HVGs don't have these genes
print("Skipping cell cycle analysis (use full gene set for cell cycle)")
print("Continuing to marker gene analysis...")

# Step 3: Neuron Marker Gene Analysis
# Updated marker genes for neurons and other cell types
marker_genes = {
    "Excitatory Neurons": ["Neurod6", "Neurod2"],
    "Glutaminergic Early Neuroblast": ["Slc17a6"],
    "Dorsal/Ventral Identity": ["Emx1", "Emx2"],
    "Endothelial Cells": ["Igfbp7", "Col4a1"],
    "GABAergic Immature Cells": ["Dlx1", "Arx"],
    "GABAergic Neuroblast": ["Gad2", "Gad1"],
    "MGE Neuroblast": ["Maf", "Lhx6", "Lhx8"],
    "MGE Progenitors": ["Nkx2-1", "Mki67"],
    "CGE Progenitors": ["Nr2f2", "Ptprz1"],
    "CGE Neuroblast": ["Prox1", "Nfib", "Dlx1", "Nr2f1", "Nrp1"],
    "LGE Neuroblast": ["Isl1", "Ebf1", "Zfhx3"]
}


# Subset marker genes to those present in the dataset
marker_genes_in_data = {}
for cell_type, markers in marker_genes.items():
    markers_found = [marker for marker in markers if marker in adata.var.index]
    marker_genes_in_data[cell_type] = markers_found


# Plot expression of marker genes for each cell type
for cell_type, markers_found in marker_genes_in_data.items():
    if markers_found: #Only proceed if there are valid markers for this cell type
        print(f"{cell_type.upper()}:")
        try:
            sc.pl.umap(
                adata,
                color=markers_found,
                #layer=layer_for_plot,
                layer=None, 
                vmin=0,
                vmax="p99",
                sort_order=False,
                frameon=False,
                cmap="inferno",
                title=f"{cell_type} res_0_25",
            )
        except ValueError as e:
            print(f"ValueError encountered: {e}. Trying alternative plotting...")
            sc.pl.umap(
                adata,
                color=[f"X_{gene}" if gene in adata.obs.columns else gene for gene in markers_found],
                vmin=0,
                vmax="p99",
                sort_order=False,
                frameon=False,
                cmap="inferno",
            )


        print("\n\n\n")
    else:
        print(f"No valid markers found for {cell_type}. Skipping...\n")


# FINAL VISUALIZATION: leiden_res0_25 + key markers
key_markers = []
for markers in marker_genes_in_data.values():
    key_markers.extend(markers[:2])  # Top 2 markers per category
key_markers = list(set(key_markers))  # Remove duplicates

print(f"\nFinal plot: leiden_res0_25 + {len(key_markers)} key markers")
sc.pl.umap(adata, color=["leiden_res0_25"] + key_markers, ncols=3, 
           vmin=0, vmax="p99", frameon=False, cmap="inferno")

# Save with leiden_res0_25
output_file = "lognorm_annotated.h5ad"
adata.write(output_file)
print(f"\n Saved: {output_file}")
print("Leiden res0_25 clusters created and ready for annotation!")