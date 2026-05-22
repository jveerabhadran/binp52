import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt

file_paths = [
    "96_6_M_KO_qc_norm.h5ad",
    "100_3_F_WT_qc_norm.h5ad",
    "114_7_F_KO_qc_norm.h5ad",
    "124_2_M_WT_qc_norm.h5ad",
    "124_4_M_KO_qc_norm.h5ad",
    "124_5_F_WT_qc_norm.h5ad",
    "148_3_M_WT_qc_norm.h5ad",
    "148_6_F_KO_qc_norm.h5ad"
]

def parse_sample_info(filename):
    parts = os.path.basename(filename).split('_')
    return {
        'sample': os.path.basename(filename).replace('_qc_norm.h5ad', ''),
        'sex': parts[2],
        'condition': parts[3]
    }

adatas = []
for file in file_paths:
    info = parse_sample_info(file)
    adata = sc.read(file)
    adata.obs['batch'] = info['sample']
    adata.obs['sex'] = info['sex']
    adata.obs['condition'] = info['condition']
    adatas.append(adata)
    print(f"Loaded {file}: {adata.n_obs} cells, {adata.n_vars} genes, {info}")

# Concatenate
adata_combined = sc.concat(
    adatas,
    axis=0,
    label="batch",
    keys=[parse_sample_info(f)['sample'] for f in file_paths],
    index_unique="-",
    join="outer"  # keep all genes from all samples
)
print(f"Combined AnnData object has {adata_combined.n_obs} cells and {adata_combined.n_vars} genes.")

all_genes = set()
total_sum = 0
for adata in adatas:
    total_sum += adata.n_vars
    all_genes.update(adata.var_names)

print(f"Sum of genes in individual files: {total_sum}")
print(f"Unique genes across all samples: {len(all_genes)}")
print(f"Genes in combined AnnData: {adata_combined.n_vars}")


# Convert sparse matrix to dense if needed (may require much memory)
if hasattr(adata_combined.X, "toarray"):
    adata_combined.X = adata_combined.X.toarray()

# Replace NaNs with zeros
adata_combined.X = np.nan_to_num(adata_combined.X, nan=0)

# Visualize batch effect before correction using PCA
sc.pp.pca(adata_combined)
sc.pl.pca_scatter(adata_combined, color='batch', title='Before batch correction')

# Batch correction with ComBat
print("Performing batch correction with ComBat...")
sc.pp.combat(adata_combined, key='batch')
print("Batch correction complete.")

# Visualize batch effect after correction
sc.pp.pca(adata_combined)
sc.pl.pca_scatter(adata_combined, color='batch', title='After batch correction')

# Save corrected object
adata_combined.write("combined_8_samples_batch_corrected.h5ad")
print("Batch-corrected combined sample file written: combined_8_samples_batch_corrected.h5ad")


