#!/usr/bin/python
"""This module plots marker genes."""

import anndata
import csv
import json
import numpy as np
import os
import pandas as pd
import scanpy as sc

from matplotlib import pyplot as plt
from pathlib import PurePath

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))


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


def parse_dge_csv(dge_csv_path) -> set[str]:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        genes = []
        for row in dge_reader:
            if float(row["FDR"]) <= 0.05:
                genes.append(row.pop(""))

        return set(genes)


print("Searching for data...")
directory_paths = set()
dge_tree = {}
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold() == "differential_gene_expression.csv":
            file_path_dge = os.path.join(root, file)
            cluster_number_path = PurePath(root).parent
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


def fill_median_table(gene, adata, cluster_key, median_table):
    clusters = adata.obs[adata_cluster_key].cat.categories
    for cluster in clusters:
        median_table.loc[gene, cluster] = np.median(
            adata[adata.obs[adata_cluster_key] == cluster].obs_vector(
                gene, layer="log1p_norm"
            )
        )


for directory_path in directory_paths:
    print(f"Processing directory {directory_path}...", flush=True)
    adata_path = os.path.join(directory_path, "cluster_relation.h5ad")
    print(f"\tLoading {adata_path}...", flush=True)
    adata = anndata.read_h5ad(adata_path)
    # Filters and sorts keys (cluster relation splitoffs) that belong to the currently processed directory.
    dge_set = set()
    keys = [x for x in dge_tree.keys() if x.startswith(directory_path)]
    keys.sort(key=sort_key_conversion)
    print(f"\tFound cluster relation tree data: {keys}", flush=True)
    for key in keys:
        print(f"\t\tProcessing {key}...", flush=True)
        output_folder_path = os.path.join(
            MOUNT_PATHS["output"],
            os.path.normpath(os.path.relpath(key, INPUT_FOLDER)),
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
        # Calculate median expression for each gene and cluster.
        adata_cluster_key = path_key_to_adata_key(key)
        clusters = adata.obs[adata_cluster_key].cat.categories

        median_dataframe = pd.DataFrame(
            np.empty((len(sorted_gene_set), len(clusters)), dtype=pd.Float64Dtype),
            columns=clusters,
            index=sorted_gene_set,
        )
        for cluster in clusters:
            print(f"\t\t\tProcessing cluster {cluster}...", flush=True)
            cluster_cells = adata[
                adata.obs[adata_cluster_key] == cluster, sorted_gene_set
            ].to_df(layer="log1p_norm")
            median_dataframe[cluster] = cluster_cells.median()

        median_dataframe.to_csv(
            os.path.join(
                output_folder_path, "differentially_expressed_genes_medians.csv"
            ),
            sep=",",
            encoding="utf-8",
        )

        # The median test has low test power, but is still preferable here
        # as only the difference in median expression not the actual distribution
        # should be tested for ordering of the clusters.