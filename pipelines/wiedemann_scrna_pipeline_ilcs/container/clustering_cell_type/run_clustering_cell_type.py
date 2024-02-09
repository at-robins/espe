#!/usr/bin/python
"""This module performs data clustering."""

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

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"] + "/"
CLUSTERING_INFO_FILE = os.path.join(MOUNT_PATHS["input"], "sample_clustering.csv")
MINIMUM_CELL_NUMBER = 20

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


def cluster_pool(sample_id: str, sample_pool: [str]):
    """
    Performs clustering.
    """
    input_files = list(
        map(
            lambda sample: os.path.join(
                INPUT_FOLDER,
                pathvalidate.sanitize_filename(sample),
                "filtered_feature_bc_matrix.h5ad",
            ),
            sample_pool,
        )
    )
    output_folder_path = os.path.join(
        MOUNT_PATHS["output"], pathvalidate.sanitize_filename(sample_id)
    )
    os.makedirs(output_folder_path, exist_ok=True)

    print(f"Processing files {input_files} as pool {sample_id}", flush=True)
    print("\tReading data...")
    adatas = list(
        map(
            anndata.read_h5ad,
            input_files,
        )
    )

    # Combines highly variant genes
    combined_deviant_genes = None
    for adata in adatas:
        if combined_deviant_genes is None:
            combined_deviant_genes = adata.var["highly_deviant"]
        else:
            combined_deviant_genes = combined_deviant_genes.combine(
                adata.var["highly_deviant"], lambda a, b: a or b, False
            )

    print("\tMerging sample pool...")
    adata_pool = anndata.concat(
        adatas,
        axis=0,
        join="outer",
        merge="same",
        uns_merge="same",
        label=None,
        keys=None,
    )
    adata_pool.var["highly_variable"] = combined_deviant_genes
    adata_pool.layers["counts"] = adata_pool.X
    np.nan_to_num(adata_pool.layers["log1p_norm"], copy=False, nan=0.0)
    adata_pool.X = adata_pool.layers["log1p_norm"]
    adata_pool.obs[CLUSTERING_GROUP_KEY] = pd.Categorical(
        np.repeat(sample_id, adata_pool.n_obs)
    )

    print(f"\tPerforming clustering for {adata_pool.n_obs} cells...", flush=True)
    if adata_pool.n_obs < MINIMUM_CELL_NUMBER:
        print("\tNot enough cells. Skipping sample...")
    else:
        sc.pp.pca(adata_pool, svd_solver="arpack", use_highly_variable=True)
        n_pcs_max = adata_pool.obsm["X_pca"].shape[1]
        if n_pcs_max < 30:
            n_pcs = n_pcs_max
        else:
            n_pcs = 30
        sc.pp.neighbors(adata_pool, n_pcs=n_pcs)
        sc.tl.umap(adata_pool)
        # output_cluster_data(
        #     adata_pool,
        #     f"{sample_id}",
        #     output_folder_path,
        # )
        sc.tl.leiden(adata_pool, key_added=CLUSTER_KEY, resolution=1.0)

        print("\tPlotting data...")
        fig = sc.pl.umap(
            adata_pool,
            color=[
                CLUSTER_KEY,
                REPLICATE_KEY,
                SAMPLE_TYPE_KEY,
            ],
            legend_loc="on data",
            show=False,
            return_fig=True,
        )
        fig.savefig(f"{output_folder_path}/umap.svg")

        print("\tWriting data to file...")
        adata_pool.write(
            f"{output_folder_path}/clustered.h5ad",
            compression="gzip",
        )

        print("\tPlotting marker genes...")
        sc.tl.rank_genes_groups(
            adata_pool,
            groupby=CLUSTER_KEY,
            method="wilcoxon",
            key_added="marker_genes_leiden",
        )

        fig = sc.pl.rank_genes_groups_dotplot(
            adata_pool,
            groupby=CLUSTER_KEY,
            standard_scale="var",
            n_genes=5,
            key="marker_genes_leiden",
            show=False,
            return_fig=True,
        )
        fig.savefig(f"{output_folder_path}/marker_genes.svg")


RESOLUTIONS_PER_ITERATION = 100
STARTING_RESOLUTION_MIN = 0.0
STARTING_RESOLUTION_MAX = 2.0


def output_cluster_data(adata: anndata.AnnData, sample_name: str, output_folder_path):
    """
    Performs clustering at different resolutions and outputs the resulting data.
    """
    print("Optimising resolution...")
    min_resolution = STARTING_RESOLUTION_MIN
    max_resolution = STARTING_RESOLUTION_MAX
    # Copies AnnData here once instead of in the actual loop where
    # needed to mitigate an associated memory leak.
    adata_tmp = adata.copy()
    cluster_data_path = os.path.join(
        output_folder_path, f"{pathvalidate.sanitize_filename(sample_name)}.csv"
    )
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
            print(f"\tTesting resolution {resolution}")
            sc.tl.leiden(adata_tmp, key_added="leiden_tmp", resolution=resolution)
            metrics_writer.writerow([resolution, *adata_tmp.obs["leiden_tmp"].values])


print("Parsing information for sample clustering...")
sample_pools = {}
with open(CLUSTERING_INFO_FILE, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.reader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        pool_id = row[0]
        pool = list(
            filter(lambda value: value != "", map(lambda value: value.strip(), row[1:]))
        )
        sample_pools[pool_id] = pool

# Clusters sample pools.
for sample_id, sample_pool in sample_pools.items():
    cluster_pool(sample_id, sample_pool)
