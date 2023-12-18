#!/usr/bin/python
"""This module performs data clustering."""

import anndata
import csv
import json
import logging
import numpy as np
import os
import pathvalidate
import scanpy as sc
import seaborn as sns

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["ilc_composition"] + "/"
CLUSTERING_INFO_FILE = os.path.join(MOUNT_PATHS["input"], "sample_clustering.csv")
CELL_TYPE_KEY = "cell_type"
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
    Performs clustering.
    """
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...")
    adata = anndata.read_h5ad(file_path_input)
    adata.layers["counts"] = adata.X
    adata.X = adata.layers["log1p_norm"]
    # Performing HGV and PCA first to reduce dimensionality for UMAP.
    adata.var["highly_variable"] = adata.var["highly_deviant"]

    cell_types = adata.obs[CELL_TYPE_KEY].cat.categories

    for cell_type in cell_types:
        print(f"\tPerform clustering for cell type {cell_type}...", flush=True)
        adata_subset = adata[adata.obs[CELL_TYPE_KEY] == cell_type].copy()
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
            sc.tl.leiden(adata_subset, key_added="leiden_res0_25", resolution=0.25)
            sc.tl.leiden(adata_subset, key_added="leiden_res0_50", resolution=0.5)
            sc.tl.leiden(adata_subset, key_added="leiden_res1_00", resolution=1.0)

            print("\tPlotting data...")
            fig = sc.pl.umap(
                adata_subset,
                color=["leiden_res0_25", "leiden_res0_50", "leiden_res1_00", "batch"],
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


def cluster_pool(sample_id: str, sample_pool: [str], metrics_writer):
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
                adata_subset, f"{sample_id}_{cell_type}", metrics_writer
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


RESOLUTIONS_PER_ITERATION = 64
STARTING_RESOLUTION_MIN = 0.1
STARTING_RESOLUTION_MAX = 2.0
MAX_ITERATIONS = 1
STABILITY_SLIDING_WINDOW = 3


def determine_optimal_clusters(
    adata: anndata.AnnData, sample_name: str, metrics_writer
):
    """
    Itertively determines the optimal resolution for clustering.
    """
    print("Optimising resolution...")
    current_iteration = 1
    min_resolution = STARTING_RESOLUTION_MIN
    max_resolution = STARTING_RESOLUTION_MAX
    # Copies AnnData here once instead of in the actual loop where
    # needed to mitigate an associated memory leak.
    adata_tmp = adata.copy()
    while current_iteration <= MAX_ITERATIONS:
        print(f"Iteration {current_iteration}/{MAX_ITERATIONS}")
        current_iteration += 1
        step_increase = (max_resolution - min_resolution) / (
            RESOLUTIONS_PER_ITERATION - 1
        )
        resolutions = [
            min_resolution + step_increase * i for i in range(RESOLUTIONS_PER_ITERATION)
        ]
        last_aggregate = None
        stabilities = []
        for resolution in resolutions:
            print(f"Testing resolution {resolution}")
            sc.tl.leiden(adata_tmp, key_added="leiden_tmp", resolution=resolution)
            # On population basis.
            current_aggregate = aggregate_clusters(adata_tmp)
            if last_aggregate is not None:
                stability = calc_stability(last_aggregate, current_aggregate)
                stabilities.append(stability)
                metrics_writer.writerow(
                    [sample_name, "Stability", resolution, stability]
                )
                print(f"\tOverall stability: {stability}")
            last_aggregate = current_aggregate
    sliding_window_normalisation = []
    for i in range(len(stabilities) - STABILITY_SLIDING_WINDOW + 1):
        current_window = stabilities[i : i + STABILITY_SLIDING_WINDOW]
        sliding_window_normalisation.append(
            sum(current_window) / STABILITY_SLIDING_WINDOW
        )
    print(f"\tNormalised stability: {sliding_window_normalisation}")


def aggregate_clusters(adata: anndata.AnnData):
    """
    Aggregates cell clusters.
    """
    print("\tAggregating cell clusters...")
    aggregate = {}
    for barcode, cluster in adata.obs["leiden_tmp"].items():
        if cluster in aggregate:
            aggregate[cluster].add(barcode)
        else:
            aggregate[cluster] = {barcode}
    return aggregate


def calc_stability(low_resolution, high_resolution):
    """
    Calculates the overall cluster stability between two resolutions.
    """
    print("\tCalulating stability...")
    cluster_stabilities = []
    for query_set in high_resolution.values():
        tmp_cluster_stability = []
        for reference_set in low_resolution.values():
            overlap = query_set.intersection(reference_set)
            tmp_cluster_stability.append(len(overlap) / len(query_set))
        # Use the maximum overlap found to calculate the cluster stability.
        cluster_stabilities.append(max(tmp_cluster_stability))
    print(
        f"\tNumber of clusters: {len(cluster_stabilities)} ; stabilities: {cluster_stabilities}",
        flush=True,
    )
    return sum(cluster_stabilities) / len(cluster_stabilities)


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
    metrics_writer.writerow(
        [
            "Sample",
            "Analysis",
            "Resolution",
            "Stability",
        ]
    )
    # Clusters sample pools.
    for sample_id, sample_pool in sample_pools.items():
        cluster_pool(sample_id, sample_pool, metrics_writer)


# Iterates over all sample directories and processes them conserving the directory structure.
# for root, dirs, files in os.walk(INPUT_FOLDER):
#     for file in files:
#         if file.casefold().endswith("filtered_feature_bc_matrix.h5ad"):
#             file_path_input = os.path.join(root, file)
#             output_folder_path = os.path.join(
#                 MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
#             )
#             os.makedirs(output_folder_path, exist_ok=True)
#             process_data(file_path_input, output_folder_path)
