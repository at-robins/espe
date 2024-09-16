#!/usr/bin/python
"""This module samples clusters at different resolutions."""

import anndata
import csv
import json
import os
import scanpy as sc


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"

RESOLUTIONS_PER_ITERATION = 201
STARTING_RESOLUTION_MIN = 0.0
STARTING_RESOLUTION_MAX = 4.0

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
    Performs clustering at different resolutions and outputs the resulting data.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)

    print("\tCalculating clustering at different resolutions...")
    min_resolution = STARTING_RESOLUTION_MIN
    max_resolution = STARTING_RESOLUTION_MAX
    # Copies AnnData here once instead of in the actual loop where
    # needed to mitigate an associated memory leak.
    adata_tmp = adata.copy()
    cluster_data_path = os.path.join(output_folder_path, "cluster_resolution_data.csv")
    with open(cluster_data_path, mode="w", newline="", encoding="utf-8") as csvfile:
        metrics_writer = csv.writer(
            csvfile,
            dialect="unix",
            delimiter=",",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
        )

        step_increase = (max_resolution - min_resolution) / (
            RESOLUTIONS_PER_ITERATION - 1
        )
        resolutions = [
            min_resolution + step_increase * i for i in range(RESOLUTIONS_PER_ITERATION)
        ]
        for resolution in resolutions:
            print(f"\tTesting resolution {resolution}...", flush=True)
            sc.tl.leiden(adata_tmp, key_added="leiden_tmp", resolution=resolution)
            metrics_writer.writerow([resolution, *adata_tmp.obs["leiden_tmp"].values])


print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("clustered.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            os.makedirs(output_folder_path, exist_ok=True)
            process_data(file_path_input, output_folder_path)
