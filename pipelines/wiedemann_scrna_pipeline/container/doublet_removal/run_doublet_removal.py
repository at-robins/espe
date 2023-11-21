#!/usr/bin/python
"""This module removes doublets."""

import anndata
import json
import os


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["doublet_detection"] + "/"


def process_data(file_path_filtered, output_folder_path):
    """
    Removes marked doublets.
    """
    print(f"Processing file {file_path_filtered}", flush=True)
    print("\tReading filtered data...")
    adata = anndata.read_h5ad(file_path_filtered)

    print("\tRemoving doublets...")
    n_cells_prefilter = adata.n_obs
    adata = adata[adata.obs["doublet_class"] == "singlet"].copy()
    n_cells_postfilter = adata.n_obs
    print(f"\tRemoved {n_cells_prefilter - n_cells_postfilter} doublets.")

    print("\tWriting filtered data to file...")
    adata.write(
        f"{output_folder_path}/filtered_feature_bc_matrix.h5ad", compression="gzip"
    )


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
            file_path_filtered = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_filtered, output_folder_path)
