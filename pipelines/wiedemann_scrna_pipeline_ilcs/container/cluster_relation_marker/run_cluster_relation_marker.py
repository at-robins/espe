#!/usr/bin/python
"""This module plots marker genes."""

import anndata
import csv
import json
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns

from matplotlib import pyplot as plt
from pathlib import PurePath

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_DGE = MOUNT_PATHS["dependencies"]["cluster_relation_dge"]
INPUT_FOLDER_DGE_MERGE = MOUNT_PATHS["dependencies"]["cluster_relation_dge_merge"]

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


def parse_dge_csv(dge_csv_path) -> set[str]:
    """
    Parses a CSV file containing differential gene expression data.
    """
    dge_file = pd.read_csv(
        dge_csv_path,
        sep=",",
        header=0,
        index_col=0,
        encoding="utf-8",
    )
    return set(dge_file.head(3).index.to_list())


print("Searching for data...")
directory_paths = set()
dge_tree = {}
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_DGE_MERGE):
    for file in files:
        if file.endswith("_merged_dge_filtered.csv"):
            file_path_dge = os.path.join(root, file)
            relative_root = os.path.normpath(
                os.path.relpath(root, INPUT_FOLDER_DGE_MERGE)
            )
            cluster_number_path = PurePath(relative_root)
            super_directory_path = cluster_number_path.parent
            directory_paths.add(str(super_directory_path))
            print(
                f"\tCollecting genes from file {file_path_dge} with {cluster_number_path.name} clusters in subdirectory {super_directory_path}...",
                flush=True,
            )
            genes = parse_dge_csv(file_path_dge)
            tree_genes_key = str(cluster_number_path)
            tree_genes = dge_tree.get(tree_genes_key)
            if tree_genes is None:
                dge_tree[tree_genes_key] = genes
            else:
                dge_tree[tree_genes_key] = tree_genes | genes


def sort_key_conversion(gene_key) -> int:
    """
    Converts the path key to the cluster number.
    """
    return int(PurePath(gene_key).name)


def path_key_to_adata_key(path_key) -> int:
    """
    Converts the path key to the respective adata key.
    """
    return f"number_of_clusters_{PurePath(path_key).name}"


for directory_path in directory_paths:
    print(f"Processing directory {directory_path}...", flush=True)
    adata_path = os.path.join(INPUT_FOLDER_DGE, directory_path, "cluster_relation.h5ad")
    print(f"\tLoading {adata_path}...", flush=True)
    adata = anndata.read_h5ad(adata_path)
    # Filters and sorts keys (cluster relation splitoffs) that
    # belong to the currently processed directory.
    dge_set = set()
    keys = [x for x in dge_tree.keys() if x.startswith(directory_path)]
    keys.sort(key=sort_key_conversion)
    print(f"\tFound cluster relation tree data: {keys}", flush=True)
    for key in keys:
        print(f"\t\tProcessing {key}...", flush=True)
        output_folder_path = os.path.join(
            MOUNT_PATHS["output"],
            os.path.normpath(key),
        )
        os.makedirs(output_folder_path, exist_ok=True)
        dge_set = dge_set | dge_tree[key]
        print("\t\tWriting total differentially expressed genes to file...", flush=True)
        sorted_gene_set = sorted(dge_set)
        with open(
            os.path.join(
                output_folder_path, "differentially_expressed_genes_total.json"
            ),
            mode="wt",
            encoding="utf-8",
        ) as gene_file:
            json.dump(sorted_gene_set, gene_file)
        print("\t\tCalculating gene medians...", flush=True)
        adata_cluster_key = path_key_to_adata_key(key)
        clusters = adata.obs[adata_cluster_key].cat.categories

        median_dataframe = pd.DataFrame(
            np.empty((len(sorted_gene_set), len(clusters)), dtype=pd.Float64Dtype),
            columns=clusters,
            index=sorted_gene_set,
        )
        cluster_dictionary = {}
        for cluster in clusters:
            print(f"\t\t\tProcessing cluster {cluster}...", flush=True)
            cluster_cells = adata[
                adata.obs[adata_cluster_key] == cluster, sorted_gene_set
            ].to_df(layer="log1p_norm")
            cluster_dictionary[cluster] = cluster_cells
            median_dataframe[cluster] = cluster_cells.mean()

        median_dataframe.to_csv(
            os.path.join(
                output_folder_path, "differentially_expressed_genes_means.csv"
            ),
            sep=",",
            encoding="utf-8",
        )

        plot = sns.clustermap(
            median_dataframe,
            cbar_kws={"label": "scaled normalised mean counts"},
            cmap="vlag",
            standard_scale=0,
            center=0.5,
            figsize=(
                6.0 + 0.6 * len(median_dataframe.columns),
                4.0 + 0.3 * len(median_dataframe.index),
            ),
        )
        plot.ax_row_dendrogram.set_visible(False)
        plot.ax_heatmap.grid(False)
        plot.savefig(
            os.path.join(
                output_folder_path,
                "clustermap_normalised.svg",
            )
        )
        reordered_dataframe = median_dataframe.iloc[
            plot.dendrogram_row.reordered_ind, plot.dendrogram_col.reordered_ind
        ]
        reordered_dataframe.to_csv(
            os.path.join(
                output_folder_path,
                "differentially_expressed_genes_normalised_means_reordered.csv",
            ),
            sep=",",
            encoding="utf-8",
        )
        plt.close()

        fig = sc.pl.stacked_violin(
            adata,
            sorted_gene_set,
            groupby=adata_cluster_key,
            swap_axes=False,
            dendrogram=True,
            show=False,
            return_fig=True,
        )
        fig.savefig(
            os.path.join(output_folder_path, "violine.svg"),
            format="svg",
        )
        plt.close()

        fig = sc.pl.dotplot(
            adata,
            sorted_gene_set,
            groupby=adata_cluster_key,
            dendrogram=True,
            show=False,
            return_fig=True,
        )
        fig.savefig(
            os.path.join(output_folder_path, "dotplot.svg"),
            format="svg",
        )
        plt.close()

        # The replicates.
        REPLICATE_KEY = "replicate_name"
        # The samples that consists of different replicates.
        SAMPLE_TYPE_KEY = "sample_type"
        # The clustering information.
        CLUSTER_KEY = "leiden_clustering"

        print("\tPlotting data...")
        all_umap_samples = [
            adata_cluster_key,
            REPLICATE_KEY,
            SAMPLE_TYPE_KEY,
            *sorted_gene_set,
        ]
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
