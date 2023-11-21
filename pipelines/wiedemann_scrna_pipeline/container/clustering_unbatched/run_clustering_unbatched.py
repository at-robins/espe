#!/usr/bin/python
"""This module performs data clustering for QC."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["feature_selection"] + "/"
UNBATCHED_SUBFOLDER = "unbatched"
INPUT_FOLDER_UNBATCHED = f"{INPUT_FOLDER}{UNBATCHED_SUBFOLDER}/"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]

def process_data(file_path_input, output_folder_path):
    """
    Performs clustering.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading filtered data...")
    adata = anndata.read_h5ad(file_path_input)

    print("\tPerforming clustering...")
    adata.X = adata.layers["log1p_norm"]
    # setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    print("\tPlotting data...")
    fig = sc.pl.umap(
        adata,
        color=["total_counts", "pct_counts_mt", "doublet_score", "doublet_class"],
        show=False,
        return_fig=True,
    )
    fig.tight_layout()
    fig.savefig(f"{output_folder_path}/umap_qc.svg")



# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_UNBATCHED):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER_UNBATCHED)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
