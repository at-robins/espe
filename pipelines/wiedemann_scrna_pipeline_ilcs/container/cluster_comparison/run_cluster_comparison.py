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

# import pertpy as pt
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

# The replicates.
REPLICATE_KEY = "replicate_name"
# The samples that consists of different replicates.
SAMPLE_TYPE_KEY = "sample_type"
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
sc.settings.figdir = ""


def is_cluster_obs(obs_name):
    """
    Returns True if the observation name contains cluster information.
    """
    return obs_name.startswith("number_of_clusters_")


def cluster_name_to_cluster(obs_name):
    """
    Returns the cluster number from the observation name.
    """
    return obs_name.removeprefix("number_of_clusters_")


def process_data(file_path_input, output_folder_path):
    """
    Computes cluster statistics.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    for obs_name in adata.obs:
        if is_cluster_obs(obs_name=obs_name):
            cluster = cluster_name_to_cluster(obs_name=obs_name)
            print(f"\tProcessing {cluster} clusters...", flush=True)
            final_output_folder = os.path.join(output_folder_path, cluster)
            os.makedirs(final_output_folder, exist_ok=True)
            cluster_statistics(
                adata=adata,
                cluster_key=obs_name,
                output_folder_path=final_output_folder,
            )


def cluster_statistics(adata, cluster_key, output_folder_path):
    """
    Computes cluster statistics.
    """
    # adata.X = adata.layers["log1p_norm"]
    # milo = pt.tl.Milo()
    # mdata = milo.load(adata)
    # milo.make_nhoods(mdata, prop=0.1)
    # nhood_size = adata.obsm["nhoods"].toarray().sum(0)
    # fig, ax = plt.subplots()
    # ax.hist(nhood_size, bins=20)
    # ax.set(xlabel="Cells in neighbourhood", ylabel="Neighbourhoods")
    # barplot_name = "neighbourhood_histogram.svg"
    # fig.savefig(
    #     os.path.join(output_folder_path, pathvalidate.sanitize_filename(barplot_name))
    # )
    # milo.count_nhoods(mdata, sample_col=REPLICATE_KEY)
    # mean_n_cells = mdata["milo"].X.toarray().mean(0)
    # fig, ax = plt.subplots()
    # ax.plot(nhood_size, mean_n_cells, ".")
    # ax.set(xlabel="Cells in neighbourhood", ylabel="Mean cells per sample in neighbourhood")
    # b_name = "neighbourhood_mean.svg"
    # fig.savefig(
    #     os.path.join(output_folder_path, pathvalidate.sanitize_filename(b_name))
    # )
    batches = adata.obs[REPLICATE_KEY].cat.categories
    samples = adata.obs[SAMPLE_TYPE_KEY].cat.categories
    clusters = adata.obs[cluster_key].cat.categories

    stat_data = {"batch": [], "sample": [], "cluster": [], "value": []}
    norm_factors = {}
    for cluster in clusters:
        print(f"\tCalculating statistics for cluster {cluster}...", flush=True)
        for batch in batches:
            adata_batch_subset = adata[adata.obs[REPLICATE_KEY] == batch]
            adata_cluster_subset = adata_batch_subset[
                adata_batch_subset.obs[cluster_key] == cluster
            ]
            stat_data["batch"].append(batch)
            stat_data["sample"].append(
                adata_batch_subset.obs[SAMPLE_TYPE_KEY].cat.categories[0]
            )
            stat_data["cluster"].append(cluster)
            frequency = adata_cluster_subset.n_obs / adata_batch_subset.n_obs
            stat_data["value"].append(frequency)
            if norm_factors.get(cluster) is None or norm_factors[cluster] < frequency:
                norm_factors[cluster] = frequency

    stat_data_frame = pd.DataFrame(data=stat_data)

    # Calculates normalisation factors.
    def normalise_values(stat_row):
        """
        Normalises the values
        """
        return stat_row["value"] / norm_factors[stat_row["cluster"]]

    stat_data_frame["norm"] = stat_data_frame.apply(normalise_values, axis=1)

    print("\tPlotting data...")
    ordering = adata.obs[SAMPLE_TYPE_KEY].cat.categories.to_numpy()
    sorted(ordering, key=str.casefold)
    barplot_name = "relative_cluster_frequency.svg"
    fig, ax = plt.subplots(figsize=(stat_data_frame.shape[0] / 4, 8))
    sns.barplot(
        stat_data_frame,
        x="cluster",
        y="norm",
        hue="sample",
        ax=ax,
        errorbar="se",
        capsize=0.15,
        hue_order=ordering,
    )
    ax.set_ylim(0, None)
    ax.set(xlabel="Cluster", ylabel="Normalised cell frequency")
    legend = ax.get_legend()
    legend.set_title("Samples")
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_folder_path, pathvalidate.sanitize_filename(barplot_name))
    )
    plt.close(fig)


# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("cluster_relation.h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            process_data(file_path_input, output_folder_path)
