#!/usr/bin/python
"""This module aggregates single cell data into pseudo-bulk data."""

import anndata
import json
import numpy as np
import os
import pandas as pd
import scanpy as sc

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustering_cell_type"] + "/"

# The replicates.
REPLICATE_KEY = "replicate_name"
# The samples that consists of different replicates.
SAMPLE_TYPE_KEY = "sample_type"
# The clustering information.
CLUSTER_KEY = "leiden_clustering"
# The sample groups that have been clustered together.
CLUSTERING_GROUP_KEY = "clustering_group"

OUTPUT_CLUSTER_KEY = "cluster"
OUTPUT_SAMPLE_KEY = "sample"
OUTPUT_CLUSTERING_GROUP_KEY = "clusteringgroup"
OUTPUT_REPLICATE_KEY = "replicate"
OUTPUT_CELL_NUMBER = "cellcount"

# Setup of scanpy.
sc.settings.verbosity = 2


def aggregate_pseudobulk_data(file_path_sample, cluster_exclude=False):
    """
    Aggregates pseudobulk data.
    """
    print(f"Processing file {file_path_sample}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_sample)
    adata.X = adata.layers["counts"]

    clusters = adata.obs[CLUSTER_KEY].cat.categories
    replicates = adata.obs[REPLICATE_KEY].cat.categories

    pseudobulk_dataframe = pd.DataFrame()
    cluster_array = []
    sample_array = []
    clustering_group_array = []
    replicate_array = []
    n_obs_array = []

    print(f"\tNumber of total observations: {adata.n_obs}")
    for cluster in clusters:
        for replicate in replicates:
            if cluster_exclude:
                print(f"\tProcessing replicate {replicate} without cluster {cluster}")
                cluster_mask = adata.obs[CLUSTER_KEY] != cluster
            else:
                print(f"\tProcessing replicate {replicate} and cluster {cluster}")
                cluster_mask = adata.obs[CLUSTER_KEY] == cluster
            replicate_mask = adata.obs[REPLICATE_KEY] == replicate
            adata_subset = adata[np.logical_and(cluster_mask, replicate_mask)]
            n_obs_subset = adata_subset.n_obs
            print(f"\t\tNumber of subset observations: {n_obs_subset}")
            if n_obs_subset > 0:
                # The replicate identifier is already specified so there can only be one sample type and clustering group.
                clustering_group = adata_subset.obs[
                    CLUSTERING_GROUP_KEY
                ].cat.categories[0]
                sample_name = adata_subset.obs[SAMPLE_TYPE_KEY].cat.categories[0]
                aggregated_row = adata_subset.to_df().agg(np.sum)
                cluster_array.append(cluster)
                sample_array.append(sample_name)
                clustering_group_array.append(clustering_group)
                replicate_array.append(replicate)
                n_obs_array.append(n_obs_subset)
                pseudobulk_dataframe = pd.concat(
                    [pseudobulk_dataframe, aggregated_row.to_frame().T],
                    ignore_index=False,
                    join="outer",
                )
            else:
                print("\t\tNo observations found. Skipping subset...")
    return (
        pseudobulk_dataframe,
        sample_array,
        clustering_group_array,
        replicate_array,
        cluster_array,
        n_obs_array,
    )


def process_data(cluster_exclude=False):
    # Aggregated pseudo-bulk scRNA-Seq data.
    pseudobulk_data = pd.DataFrame()
    samples = []
    clustering_groups = []
    replicates = []
    clusters = []
    n_cells = []

    # Iterates over all sample directories and collects pseudo bulk data.
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.casefold().endswith("clustered.h5ad"):
                input_file_path = os.path.join(root, file)
                (
                    tmp_pseudobulk,
                    tmp_samples,
                    tmp_clustering_groups,
                    tmp_replicates,
                    tmp_clusters,
                    tmp_n_cells,
                ) = aggregate_pseudobulk_data(
                    file_path_sample=input_file_path,
                    cluster_exclude=cluster_exclude,
                )
                pseudobulk_data = pd.concat(
                    [pseudobulk_data, tmp_pseudobulk],
                    ignore_index=False,
                    join="outer",
                )
                samples.extend(tmp_samples)
                clustering_groups.extend(tmp_clustering_groups)
                replicates.extend(tmp_replicates)
                clusters.extend(tmp_clusters)
                n_cells.extend(tmp_n_cells)

    # Replaces NAs with zeros.
    pseudobulk_data.fillna(0, inplace=True)

    # Sets row names.
    row_names = list(
        map(
            lambda name: name.replace(" ", "_"),
            [
                f"{sample}_{replicate}_{cluster}"
                for sample, replicate, cluster in zip(samples, replicates, clusters)
            ],
        )
    )
    pseudobulk_data.index = row_names

    # Sets additional variables.
    pseudobulk_adata = anndata.AnnData(X=pseudobulk_data)
    pseudobulk_adata.obs[OUTPUT_SAMPLE_KEY] = samples
    pseudobulk_adata.obs[OUTPUT_CLUSTERING_GROUP_KEY] = clustering_groups
    pseudobulk_adata.obs[OUTPUT_REPLICATE_KEY] = replicates
    pseudobulk_adata.obs[OUTPUT_CLUSTER_KEY] = clusters
    pseudobulk_adata.obs[OUTPUT_CELL_NUMBER] = n_cells

    # Sets categorical data to a categorical data type.
    pseudobulk_adata.obs[OUTPUT_SAMPLE_KEY] = pseudobulk_adata.obs[
        OUTPUT_SAMPLE_KEY
    ].astype("category")
    pseudobulk_adata.obs[OUTPUT_REPLICATE_KEY] = pseudobulk_adata.obs[
        OUTPUT_REPLICATE_KEY
    ].astype("category")
    pseudobulk_adata.obs[OUTPUT_CLUSTER_KEY] = pseudobulk_adata.obs[
        OUTPUT_CLUSTER_KEY
    ].astype("category")

    print("\tWriting pseudobulk data to file...")
    output_file_name = (
        "cluster_exclude_pseudobulk.h5ad"
        if cluster_exclude
        else "cluster_pseudobulk.h5ad"
    )
    pseudobulk_adata.write(
        os.path.join(MOUNT_PATHS["output"], output_file_name), compression="gzip"
    )


process_data(cluster_exclude=False)
process_data(cluster_exclude=True)
