#!/usr/bin/python
"""This module shows marker gene expression in the clustered data."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pathvalidate
import scanpy as sc
import seaborn as sns
import sys

from matplotlib import pyplot as plt

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"
MARKER_FILE = os.path.join(MOUNT_PATHS["input"], "markers.csv")

# The replicates.
REPLICATE_KEY = "replicate_name"
# The samples that consists of different replicates.
SAMPLE_TYPE_KEY = "sample_type"
# The clustering information.
CLUSTER_KEY = "leiden_clustering"

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
    Marks the specified gene.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    print("\tReading markers...")
    genes = []
    with open(MARKER_FILE, newline="", encoding="utf-8") as csvfile:
        info_reader = csv.reader(csvfile, dialect="unix", delimiter=",", quotechar='"')
        for row in info_reader:
            gene = row[0]
            genes.append(gene)
            if gene not in adata.var_names:
                print(
                    f"\tMarker gene {gene} not present in the data. Adding with value of NaN..."
                )
                adata.obs[gene] = [np.nan] * adata.n_obs
    print(f"\tUsing  markers: {genes}")

    print("\tPlotting samples...")
    samples = adata.obs[SAMPLE_TYPE_KEY].cat.categories
    for sample in samples:
        fig = sc.pl.umap(
            adata,
            color=[SAMPLE_TYPE_KEY],
            groups=sample,
            wspace=1,
            show=False,
            return_fig=True,
        )
        fig.savefig(
            os.path.join(
                output_folder_path,
                pathvalidate.sanitize_filename(f"sample_{sample}.svg"),
            )
        )
        plt.close(fig)

    print("\tPlotting data...")
    all_umap_samples = [CLUSTER_KEY, REPLICATE_KEY, SAMPLE_TYPE_KEY, *genes]
    # Plots
    chunk_size = 4 * 8
    for chunk_index in range(0, len(all_umap_samples), chunk_size):
        fig = sc.pl.umap(
            adata,
            color=all_umap_samples[chunk_index : chunk_index + chunk_size],
            wspace=1,
            show=False,
            return_fig=True,
        )
        fig.savefig(
            os.path.join(output_folder_path, f"marker_expression_{chunk_index}.svg")
        )
        plt.close(fig)


if not os.path.isfile(MARKER_FILE):
    print("No marker file present. Exiting...")
    sys.exit()

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
