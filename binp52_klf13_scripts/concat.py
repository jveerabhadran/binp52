import numpy as np
import scanpy as sc
import scipy.sparse as sp
import os

file_paths = [
    "binp52_klf13_ht/1_1_M_HT.h5ad",
    "binp52_klf13_ht/1_2_F_HT.h5ad",
    "binp52_klf13_ht/1_3_M_HT.h5ad",
    "binp52_klf13_ht/1_4_F_HT.h5ad",
    "binp52_klf13_ht/2_1_M_WT.h5ad",
    "binp52_klf13_ht/2_2_F_WT.h5ad",
    "binp52_klf13_ht/2_3_M_WT.h5ad",
    "binp52_klf13_ht/2_4_F_WT.h5ad",
]

def parse_sample_info(filename):
    parts = os.path.basename(filename).split("_")
    return {
        "sample": os.path.basename(filename).replace(".h5ad", ""),
        "sex": parts[2],
        "condition": parts[3],
    }

def make_combined_all_layers():
    # Load first file to see all available layers
    first_file = file_paths[0]
    adata_sample = sc.read(first_file)
    all_layers = list(adata_sample.layers.keys())
    print(f"Found layers: {all_layers}")

    # Prepare adatas with all layers preserved
    adatas = []
    for file in file_paths:
        info = parse_sample_info(file)
        adata = sc.read(file)

        # Keep original X as counts (or log1p_norm), copy all layers
        adata.obs["batch"] = info["sample"]
        adata.obs["sex"] = info["sex"]
        adata.obs["condition"] = info["condition"]
        adatas.append(adata)
        print(f"Loaded {file}: {adata.n_obs} cells, {adata.n_vars} genes, {len(adata.layers)} layers")

    # Concatenate keeping all layers
    adata_combined = sc.concat(
        adatas,
        axis=0,
        label="batch",
        keys=[parse_sample_info(f)["sample"] for f in file_paths],
        index_unique="-",
        join="outer",
    )

    # Handle NaNs in X (preserve layers unchanged)
    if sp.issparse(adata_combined.X):
        adata_combined.X.data[np.isnan(adata_combined.X.data)] = 0
    else:
        adata_combined.X = np.nan_to_num(adata_combined.X, nan=0)

    # Save with all layers
    out_file = "combined_8_samples_all_layers_outer.h5ad"
    adata_combined.write(out_file)
    print(f"Saved {out_file} with {len(adata_combined.layers)} layers")
    print(f"Final shape: {adata_combined.shape}")

if __name__ == "__main__":
    make_combined_all_layers()
