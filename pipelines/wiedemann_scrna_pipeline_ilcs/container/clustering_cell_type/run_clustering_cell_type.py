#!/usr/bin/python
"""This module performs data clustering."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["ilc_composition"] + "/"
DEFAULT_BATCH_KEY = "cell_type"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    scanpy=True,
    # # In case of bitmap exports, use high quality.
    dpi=300,
    dpi_save=300,
    # Export as SVG.
    format="svg",
    vector_friendly=True,
    # Use transparent background.
    transparent=True,
    facecolor=None,
    # Remove frames.
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


def process_data(file_path_input, output_folder_path):
    """
    Performs clustering.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    adata.X = adata.layers["log1p_norm"]
    # Performing HGV and PCA first to reduce dimensionality for UMAP.
    adata.var["highly_variable"] = adata.var["highly_deviant"]

    batches = adata.obs[DEFAULT_BATCH_KEY].cat.categories

    for batch in batches:
        print(f"\tPerform clustering for cell type {batch}...", flush=True)
        adata_subset = adata[adata.obs[DEFAULT_BATCH_KEY] == batch].copy()
        if adata_subset.n_obs < 20:
            print("\tNot enough cells. Skipping cell type...")
        else:
            sc.pp.pca(adata_subset, svd_solver="arpack", use_highly_variable=True)
            n_pcs_max = adata_subset.obsm["X_pca"].shape[1]
            if n_pcs_max < 30:
                n_pcs = n_pcs_max
            else:
                n_pcs = 30
            sc.pp.neighbors(adata_subset, n_pcs=n_pcs)
            sc.tl.umap(adata_subset)
            sc.tl.leiden(adata_subset, key_added="leiden_res0_25", resolution=0.25)
            sc.tl.leiden(adata_subset, key_added="leiden_res0_5", resolution=0.5)
            sc.tl.leiden(adata_subset, key_added="leiden_res1", resolution=1.0)

            print("\tPlotting data...")
            fig = sc.pl.umap(
                adata_subset,
                color=["leiden_res0_25", "leiden_res0_5", "leiden_res1", "batch"],
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            fig.savefig(f"{output_folder_path}/umap_{batch}.svg")


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
