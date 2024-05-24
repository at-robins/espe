#!/usr/bin/python
"""This module plots marker genes."""

import anndata
import csv
import json
import numpy as np
import os
import pandas as pd
import pathvalidate
import scanpy as sc
import scipy
import seaborn as sns

from matplotlib import pyplot as plt
from pathlib import PurePath

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_DGE = MOUNT_PATHS["dependencies"]["cluster_relation_dge"] + "/"
INPUT_FOLDER_TREE = MOUNT_PATHS["dependencies"]["cluster_relation_tree"] + "/"
SIGNIFICANCE_CUTOFF = 0.05


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


def parse_dge_csv(dge_csv_path, depth) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        genes = []
        lfc = []
        significance = []
        for row in dge_reader:
            genes.append(row[""])
            lfc.append(float(row["logFC"]))
            significance.append(
                depth if float(row["FDR"]) <= SIGNIFICANCE_CUTOFF else -1
            )

        return pd.DataFrame(
            data={"logFC": lfc, "significant": significance}, index=genes
        )


def merge_dges(dge_a, dge_b):
    """
    Merges two differential expression dataframes.
    """
    return pd.DataFrame(
        data={
            "logFC": dge_a["logFC"].add(dge_b["logFC"], fill_value=0),
            "significant": dge_a["significant"].combine(
                dge_b["significant"], max, fill_value=-1
            ),
        },
        index=dge_a.index,
    )


def load_and_merge_dges(dge_paths):
    """
    Loads the specifed DGE files and merges them.
    The files need to be sorted by increasing depth.
    """
    merged_dges = parse_dge_csv(dge_csv_path=dge_paths.pop(0), depth=0)
    for depth, dge_path in enumerate(dge_paths, start=1):
        merged_dges = merge_dges(
            merged_dges, parse_dge_csv(dge_csv_path=dge_path, depth=depth)
        )
    return merged_dges


def cluster_to_merged_dges(number_of_clusters, cluster, cluster_tree, sub_folder):
    """ """
    # Filters out unnecessary information of higher order clusterings.
    filtered_tree = [
        x for x in cluster_tree if x["number_of_clusters"] < number_of_clusters
    ]
    current_cluster = cluster
    current_noc = number_of_clusters
    dge_file_paths = []
    for cluster_entry in filtered_tree:
        # Check all child clusters of the supercluster for the current cluster.
        for node in cluster_entry["nodes"]:
            if current_cluster in node["child_clusters"]:
                if len(node["child_clusters"]) > 1:
                    reference_cluster_ids = [
                        z for z in node["child_clusters"] if z != current_cluster
                    ]
                    dge_file_paths.append(
                        os.path.join(
                            INPUT_FOLDER_DGE,
                            sub_folder,
                            pathvalidate.sanitize_filename(str(current_noc)),
                            pathvalidate.sanitize_filename(
                                f"{current_cluster}__vs__{'_'.join(map(str, reference_cluster_ids))}"
                            ),
                            "differential_gene_expression.csv",
                        )
                    )
                current_cluster = node["cluster_id"]
                break
        current_noc = cluster_entry["number_of_clusters"]

    dge_file_paths.reverse()
    return dge_file_paths


def load_tree(tree_path):
    """
    Loads the cluster relation tree file.
    """
    print("\tLoading cluster relation tree...", flush=True)
    with open(tree_path, mode="rt", encoding="utf-8") as tree:
        tree = json.load(tree)
        tree.sort(reverse=True, key=lambda x: x["number_of_clusters"])
        print(tree, flush=True)
        return tree


# def load_counts(counts_path, tree, output_folder):
#     """
#     Loads the count matrix and plots the sampled clusters.
#     """
#     print("\tLoading count data...", flush=True)
#     adata = anndata.read_h5ad(counts_path)
#     entry_obs_names = []
#     for entry in tree:
#         entry_obs_name = n_cluster_obs_name(entry["number_of_clusters"])
#         adata.obs[entry_obs_name] = pd.Categorical(entry["clustering"])
#         entry_obs_names.append(entry_obs_name)

#     print("\tSaving count data...", flush=True)
#     adata.write(os.path.join(output_folder, "cluster_relation.h5ad"))

#     print("\tPlotting data...")
#     fig = sc.pl.umap(
#         adata,
#         color=entry_obs_names,
#         legend_loc="on data",
#         show=False,
#         return_fig=True,
#     )
#     fig.savefig(os.path.join(output_folder, "umap.svg"))
#     plt.close(fig)

#     adata.X = adata.layers["counts"]
#     return adata


# def n_cluster_obs_name(number_of_clusters):
#     """
#     Returns the observation name for the clustering info based on the number of clusters specified.
#     """dusp
#     return f"number_of_clusters_{number_of_clusters}"


print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_TREE):
    for file in files:
        if file.casefold().endswith("genealogy_cluster_resolution_data.json"):
            file_path_tree = os.path.join(root, file)
            relation_tree = load_tree(file_path_tree)
            relative_sub_directory = os.path.normpath(
                os.path.relpath(root, INPUT_FOLDER_TREE)
            )
            folder_path_dge = os.path.join(
                INPUT_FOLDER_DGE,
                relative_sub_directory,
            )
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                relative_sub_directory,
            )
            for cluster_entry in relation_tree:
                for node in cluster_entry["nodes"]:
                    dge_files = cluster_to_merged_dges(
                        number_of_clusters=cluster_entry["number_of_clusters"],
                        cluster=node["cluster_id"],
                        cluster_tree=relation_tree,
                        sub_folder=relative_sub_directory,
                    )
                    print(dge_files, flush=True)
                    merged_dge_file = load_and_merge_dges(dge_paths=dge_files)
                    final_output_folder = os.path.join(
                        output_folder_path,
                        pathvalidate.sanitize_filename(str(cluster_entry["number_of_clusters"])),
                    )
                    os.makedirs(final_output_folder, exist_ok=True)
                    final_output_file = os.path.join(
                        final_output_folder,
                        pathvalidate.sanitize_filename(
                            f"{node['cluster_id']}_merged_dge.csv"
                        ),
                    )
                    merged_dge_file.to_csv(
                        final_output_file,
                        sep=",",
                        encoding="utf-8",
                    )
#             print(
#                 f"Processing files {file_path_counts}, {file_path_sampling} and {file_path_tree}",
#                 flush=True,
#             )
#             relation_tree = load_tree(file_path_tree)
#             load_sampling(file_path_sampling, relation_tree)
#             counts = load_counts(
#                 counts_path=file_path_counts,
#                 tree=relation_tree,
#                 output_folder=output_folder_path,
#             )
#             dge_for_all_splits(
#                 adata=counts, tree=relation_tree, output_path=output_folder_path
#             )
