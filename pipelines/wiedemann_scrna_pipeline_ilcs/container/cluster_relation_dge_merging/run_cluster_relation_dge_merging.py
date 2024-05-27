#!/usr/bin/python
"""
This module merges differential gene expression data from
the cluster relation tree.
"""

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
KEY_ENRICHMENT_SCORE = "enrichment"
KEY_SIGNIFICANCE_DEPTH = "significance"


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
            data={KEY_ENRICHMENT_SCORE: lfc, KEY_SIGNIFICANCE_DEPTH: significance}, index=genes
        )


def merge_dges(dge_a, dge_b):
    """
    Merges two differential expression dataframes.
    """
    return pd.DataFrame(
        data={
            KEY_ENRICHMENT_SCORE: dge_a[KEY_ENRICHMENT_SCORE].add(dge_b[KEY_ENRICHMENT_SCORE], fill_value=0),
            KEY_SIGNIFICANCE_DEPTH: dge_a[KEY_SIGNIFICANCE_DEPTH].combine(
                dge_b[KEY_SIGNIFICANCE_DEPTH], max, fill_value=-1
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


def get_branching_dges(number_of_clusters, cluster, cluster_tree, sub_folder):
    """
    Calculate branches of the cluster relation tree for the specified input cluster
    and select the relevant files for merging of the expression data.
    """
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

    # Orders the splits from root to leaf.
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
        return tree


print("Searching for data...")
# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER_TREE):
    for file in files:
        if file.casefold().endswith("genealogy_cluster_resolution_data.json"):
            file_path_tree = os.path.join(root, file)
            print(f"Processing cluster relation tree data {file_path_tree}...")
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
            for cluster_entry_index, cluster_entry in enumerate(relation_tree):
                if cluster_entry_index == len(relation_tree) - 1:
                    # Skips the step at the root of the cluster relation tree
                    # as no information on previous branches is available.
                    break
                number_of_clusters = cluster_entry["number_of_clusters"]
                print(f"\tProcessing clustering with {number_of_clusters} clusters...")
                for node in cluster_entry["nodes"]:
                    cluster_id = node["cluster_id"]
                    print(f"\t\tProcessing cluster {cluster_id}...")
                    dge_files = get_branching_dges(
                        number_of_clusters=number_of_clusters,
                        cluster=cluster_id,
                        cluster_tree=relation_tree,
                        sub_folder=relative_sub_directory,
                    )
                    print(f"\t\tThe following DGE files have been selected: {dge_files}")
                    print("\t\tMerging gene expression data...")
                    merged_dge_file = load_and_merge_dges(dge_paths=dge_files)
                    final_output_folder = os.path.join(
                        output_folder_path,
                        pathvalidate.sanitize_filename(str(number_of_clusters)),
                    )
                    print("\t\tWriting output to file...")
                    os.makedirs(final_output_folder, exist_ok=True)
                    final_output_file = os.path.join(
                        final_output_folder,
                        pathvalidate.sanitize_filename(
                            f"{cluster_id}_merged_dge.csv"
                        ),
                    )
                    merged_dge_file.to_csv(
                        final_output_file,
                        sep=",",
                        encoding="utf-8",
                    )
