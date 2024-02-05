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
CELL_TYPE_KEY = "cell_type"
CLUSTER_KEY = "leiden_res1_00"
SAMPLE_KEY = "sample"
BATCH_KEY = "batch"
MINIMUM_CELL_NUMBER = 20

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
    # The data has already been separated by cell type.
    cell_type = adata.obs[CELL_TYPE_KEY].cat.categories[0]
    batches = adata.obs[BATCH_KEY].cat.categories
    samples = adata.obs[SAMPLE_KEY].cat.categories
    clusters = adata.obs[CLUSTER_KEY].cat.categories
    print(f"\tCell type: {cell_type}", flush=True)
    print(f"\tBatches: {batches}", flush=True)
    print(f"\tSamples: {samples}", flush=True)
    print(f"\tClusters: {clusters}", flush=True)

    stat_data = {"batch": [], "sample": [], "cluster": [], "value": []}
    for cluster in clusters:
        print(f"\tCalculating statistics for cluster {cluster}...", flush=True)
        for batch in batches:
            adata_batch_subset = adata[adata.obs[BATCH_KEY] == batch]
            adata_cluster_subset = adata_batch_subset[
                adata_batch_subset.obs[CLUSTER_KEY] == cluster
            ]
            stat_data["batch"].append(batch)
            stat_data["sample"].append(
                adata_batch_subset.obs[SAMPLE_KEY].cat.categories[0]
            )
            stat_data["cluster"].append(cluster)
            stat_data["value"].append(
                adata_cluster_subset.n_obs / adata_batch_subset.n_obs
            )

    stat_data_frame = pd.DataFrame(data=stat_data)
    print(f"\tStats: {stat_data_frame}", flush=True)
    print("\tPlotting data...")
    ordering = adata.obs[SAMPLE_KEY].cat.categories.to_numpy()
    sorted(ordering, key=str.casefold)
    barplot_name = f"relative_cluster_frequency_{cell_type}.svg"
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
    ax.set_title(cell_type)
    ax.set_ylim(0, None)
    ax.set(xlabel="Cluster", ylabel="Cell frequency")
    legend = ax.get_legend()
    legend.set_title("Samples")
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_folder_path, pathvalidate.sanitize_filename(barplot_name))
    )

    # adata.layers["counts"] = adata.X
    # adata.X = adata.layers["log1p_norm"]
    # # Performing HGV and PCA first to reduce dimensionality for UMAP.
    # adata.var["highly_variable"] = adata.var["highly_deviant"]

    # cell_types = adata.obs[CELL_TYPE_KEY].cat.categories

    # for cell_type in cell_types:
    #     print(f"\tPerform clustering for cell type {cell_type}...", flush=True)
    #     adata_subset = adata[adata.obs[CELL_TYPE_KEY] == cell_type].copy()
    #     if adata_subset.n_obs < MINIMUM_CELL_NUMBER:
    #         print("\tNot enough cells. Skipping cell type...")
    #     else:
    #         sc.pp.pca(adata_subset, svd_solver="arpack", use_highly_variable=True)
    #         n_pcs_max = adata_subset.obsm["X_pca"].shape[1]
    #         if n_pcs_max < 30:
    #             n_pcs = n_pcs_max
    #         else:
    #             n_pcs = 30
    #         sc.pp.neighbors(adata_subset, n_pcs=n_pcs)
    #         sc.tl.umap(adata_subset)
    #         sc.tl.leiden(adata_subset, key_added="leiden_res0_25", resolution=0.25)
    #         sc.tl.leiden(adata_subset, key_added="leiden_res0_50", resolution=0.5)
    #         sc.tl.leiden(adata_subset, key_added="leiden_res1_00", resolution=1.0)

    #         print("\tPlotting data...")
    #         fig = sc.pl.umap(
    #             adata_subset,
    #             color=["leiden_res0_25", "leiden_res0_50", "leiden_res1_00", "batch"],
    #             legend_loc="on data",
    #             show=False,
    #             return_fig=True,
    #         )
    #         sanitised_cell_type = pathvalidate.sanitize_filename(cell_type)
    #         fig.savefig(f"{output_folder_path}/umap_{sanitised_cell_type}.svg")

    #         print("\tWriting data to file...")
    #         adata_subset.write(
    #             f"{output_folder_path}/clustered_{sanitised_cell_type}.h5ad",
    #             compression="gzip",
    #         )

    #         print("\tPlotting marker genes...")
    #         sc.tl.rank_genes_groups(
    #             adata_subset,
    #             groupby="leiden_res1_00",
    #             method="wilcoxon",
    #             key_added="marker_genes_leiden_res1_00",
    #         )

    #         fig = sc.pl.rank_genes_groups_dotplot(
    #             adata_subset,
    #             groupby="leiden_res1_00",
    #             standard_scale="var",
    #             n_genes=5,
    #             key="marker_genes_leiden_res1_00",
    #             show=False,
    #             return_fig=True,
    #         )
    #         fig.savefig(
    #             f"{output_folder_path}/marker_genes_resolution_1_00_{sanitised_cell_type}.svg"
    #         )


def cluster_pool(sample_id: str, sample_pool: [str], metrics_writer):
    """
    Computes cluster statistics.
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
        label="sample",
        keys=sample_pool,
    )
    adata_pool.var["highly_variable"] = combined_deviant_genes
    adata_pool.layers["counts"] = adata_pool.X
    np.nan_to_num(adata_pool.layers["log1p_norm"], copy=False, nan=0.0)
    adata_pool.X = adata_pool.layers["log1p_norm"]

    cell_types = adata_pool.obs[CELL_TYPE_KEY].cat.categories
    for cell_type in cell_types:
        print(f"\tPerform clustering for cell type {cell_type}...", flush=True)
        print(adata_pool.n_obs)
        adata_subset = adata_pool[adata_pool.obs[CELL_TYPE_KEY] == cell_type].copy()
        print(adata_subset.n_obs)
        if adata_subset.n_obs < MINIMUM_CELL_NUMBER:
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
            determine_optimal_clusters(
                adata_subset,
                f"{sample_id}_{cell_type}",
                output_folder_path,
            )
            sc.tl.leiden(adata_subset, key_added="leiden_res0_25", resolution=0.25)
            sc.tl.leiden(adata_subset, key_added="leiden_res0_50", resolution=0.5)
            sc.tl.leiden(adata_subset, key_added="leiden_res1_00", resolution=1.0)

            print("\tPlotting data...")
            fig = sc.pl.umap(
                adata_subset,
                color=[
                    "leiden_res0_25",
                    "leiden_res0_50",
                    "leiden_res1_00",
                    "batch",
                    "sample",
                ],
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            sanitised_cell_type = pathvalidate.sanitize_filename(cell_type)
            fig.savefig(f"{output_folder_path}/umap_{sanitised_cell_type}.svg")

            print("\tWriting data to file...")
            adata_subset.write(
                f"{output_folder_path}/clustered_{sanitised_cell_type}.h5ad",
                compression="gzip",
            )

            print("\tPlotting marker genes...")
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby="leiden_res1_00",
                method="wilcoxon",
                key_added="marker_genes_leiden_res1_00",
            )

            fig = sc.pl.rank_genes_groups_dotplot(
                adata_subset,
                groupby="leiden_res1_00",
                standard_scale="var",
                n_genes=5,
                key="marker_genes_leiden_res1_00",
                show=False,
                return_fig=True,
            )
            fig.savefig(
                f"{output_folder_path}/marker_genes_resolution_1_00_{sanitised_cell_type}.svg"
            )


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

with open(
    f"{MOUNT_PATHS['output']}/metrics.csv", mode="w", newline="", encoding="utf-8"
) as csvfile:
    metrics_writer = csv.writer(
        csvfile, dialect="unix", delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
    )
    # metrics_writer.writerow(
    #     [
    #         "Sample",
    #         "Analysis",
    #         "Resolution",
    #         "Stability",
    #     ]
    # )
    # Clusters sample pools.
    for sample_id, sample_pool in sample_pools.items():
        # cluster_pool(sample_id, sample_pool, metrics_writer)
        sample_input_folder = os.path.join(
            INPUT_FOLDER, pathvalidate.sanitize_filename(sample_id)
        )
        # Iterates over all sample directories and processes them conserving the directory structure.
        for root, dirs, files in os.walk(sample_input_folder):
            for file in files:
                if file.startswith("clustered_") and file.endswith(".h5ad"):
                    file_path_input = os.path.join(root, file)
                    output_folder_path = os.path.join(
                        MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
                    )
                    os.makedirs(output_folder_path, exist_ok=True)
                    process_data(file_path_input, output_folder_path)
