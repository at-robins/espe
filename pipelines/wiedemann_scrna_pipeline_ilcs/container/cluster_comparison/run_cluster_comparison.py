#!/usr/bin/python
"""This module compares clusters of different samples."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pandas as pd
import pathvalidate
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"
CLUSTERING_INFO_FILE = os.path.join(MOUNT_PATHS["input"], "sample_clustering.csv")

# The replicates.
REPLICATE_KEY = "replicate_name"
# The samples that consists of different replicates.
SAMPLE_TYPE_KEY = "sample_type"
# The clustering information.
CLUSTER_KEY = "leiden_clustering"
# The sample groups that have been clustered together.
CLUSTERING_GROUP_KEY = "clustering_group"

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
    Computes cluster statistics.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    batches = adata.obs[REPLICATE_KEY].cat.categories
    samples = adata.obs[SAMPLE_TYPE_KEY].cat.categories
    clusters = adata.obs[CLUSTER_KEY].cat.categories

    stat_data = {"batch": [], "sample": [], "cluster": [], "value": []}
    for cluster in clusters:
        print(f"\tCalculating statistics for cluster {cluster}...", flush=True)
        for batch in batches:
            adata_batch_subset = adata[adata.obs[REPLICATE_KEY] == batch]
            adata_cluster_subset = adata_batch_subset[
                adata_batch_subset.obs[CLUSTER_KEY] == cluster
            ]
            stat_data["batch"].append(batch)
            stat_data["sample"].append(
                adata_batch_subset.obs[SAMPLE_TYPE_KEY].cat.categories[0]
            )
            stat_data["cluster"].append(cluster)
            stat_data["value"].append(
                adata_cluster_subset.n_obs / adata_batch_subset.n_obs
            )

    stat_data_frame = pd.DataFrame(data=stat_data)
    print("\tPlotting data...")
    ordering = adata.obs[SAMPLE_TYPE_KEY].cat.categories.to_numpy()
    sorted(ordering, key=str.casefold)
    barplot_name = "relative_cluster_frequency.svg"
    fig, ax = plt.subplots(figsize=(stat_data_frame.shape[0] / 4, 8))
    sns.barplot(
        stat_data_frame,
        x="cluster",
        y="value",
        hue="sample",
        ax=ax,
        errorbar="se",
        capsize=0.15,
        hue_order=ordering,
    )
    ax.set_ylim(0, None)
    ax.set(xlabel="Cluster", ylabel="Cell frequency")
    legend = ax.get_legend()
    legend.set_title("Samples")
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_folder_path, pathvalidate.sanitize_filename(barplot_name))
    )



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
