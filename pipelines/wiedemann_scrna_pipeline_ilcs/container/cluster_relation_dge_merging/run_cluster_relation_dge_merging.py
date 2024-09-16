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
KEY_ENRICHMENT_SCORE_RAW = "raw enrichment"
KEY_SIGNIFICANCE_DEPTH = "significance"
KEY_ENRICHMENT_SCORE_WEIGHTED = "weighted enrichment"


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


def parse_dge_csv(
    dge_csv_path, depth, max_depth, cluster_size_weighting
) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        genes = []
        lfc = []
        lfc_weighted = []
        significance = []
        for row in dge_reader:
            genes.append(row[""])
            raw_lfc = float(row["logFC"])
            depth_weighting = 2.0 ** (depth - max_depth)
            lfc.append(raw_lfc)
            lfc_weighted.append(depth_weighting * cluster_size_weighting * raw_lfc)
            significance.append(
                depth if float(row["FDR"]) <= SIGNIFICANCE_CUTOFF else -1
            )

        return pd.DataFrame(
            data={
                KEY_ENRICHMENT_SCORE_RAW: lfc,
                KEY_SIGNIFICANCE_DEPTH: significance,
                KEY_ENRICHMENT_SCORE_WEIGHTED: lfc_weighted,
            },
            index=genes,
        )


def merge_dges(dge_a, dge_b):
    """
    Merges two differential expression dataframes.
    """
    return pd.DataFrame(
        data={
            KEY_ENRICHMENT_SCORE_RAW: dge_a[KEY_ENRICHMENT_SCORE_RAW].add(
                dge_b[KEY_ENRICHMENT_SCORE_RAW], fill_value=0
            ),
            KEY_SIGNIFICANCE_DEPTH: dge_a[KEY_SIGNIFICANCE_DEPTH].combine(
                dge_b[KEY_SIGNIFICANCE_DEPTH], max, fill_value=-1
            ),
            KEY_ENRICHMENT_SCORE_WEIGHTED: dge_a[KEY_ENRICHMENT_SCORE_WEIGHTED].add(
                dge_b[KEY_ENRICHMENT_SCORE_WEIGHTED], fill_value=0
            ),
        },
        index=dge_a.index,
    )


def load_and_merge_dges(dge_paths, cluster_size_weights):
    """
    Loads the specifed DGE files and merges them.
    The files need to be sorted by increasing depth.
    """
    max_depth = len(dge_paths) - 1
    merged_dges = parse_dge_csv(
        dge_csv_path=dge_paths.pop(0),
        depth=0,
        max_depth=max_depth,
        cluster_size_weighting=cluster_size_weights.pop(0),
    )
    for depth, dge_path in enumerate(dge_paths, start=1):
        merged_dges = merge_dges(
            merged_dges,
            parse_dge_csv(
                dge_csv_path=dge_path,
                depth=depth,
                max_depth=max_depth,
                cluster_size_weighting=cluster_size_weights.pop(0),
            ),
        )
    return merged_dges


def get_branching_dges(number_of_clusters, cluster, cluster_tree, sub_folder):
    """
    Calculate branches of the cluster relation tree for the specified input cluster
    and select the relevant files for merging of the expression data.
    """
    # Filters out unnecessary information of higher order clusterings.
    filtered_tree = [
        x for x in cluster_tree if x["number_of_clusters"] <= number_of_clusters
    ]
    # The entry that contains size information of the child clusters.
    child_cluster_entry = filtered_tree.pop(0)
    current_cluster = cluster
    current_noc = number_of_clusters
    dge_file_paths = []
    cluster_size_weightings = []
    for cluster_entry in filtered_tree:
        # Check all child clusters of the supercluster for the current cluster.
        for node in cluster_entry["nodes"]:
            if current_cluster in node["child_clusters"]:
                if len(node["child_clusters"]) > 1:
                    reference_cluster_ids = [
                        z for z in node["child_clusters"] if z != current_cluster
                    ]
                    # Determines cluster size weighting.
                    sample_cluster_size = sum(
                        [
                            z["number_of_cells"]
                            for z in child_cluster_entry["nodes"]
                            if z["cluster_id"] == current_cluster
                        ]
                    )
                    reference_cluster_size = sum(
                        [
                            z["number_of_cells"]
                            for z in child_cluster_entry["nodes"]
                            if z["cluster_id"] in reference_cluster_ids
                        ]
                    )
                    size_weighting = reference_cluster_size / (
                        reference_cluster_size + sample_cluster_size
                    )
                    cluster_size_weightings.append(size_weighting)
                    # Determines DGE file path.
                    potential_dge_path = os.path.join(
                        INPUT_FOLDER_DGE,
                        sub_folder,
                        pathvalidate.sanitize_filename(str(current_noc)),
                        pathvalidate.sanitize_filename(
                            f"{current_cluster}__vs__{'_'.join(map(str, reference_cluster_ids))}"
                        ),
                        "differential_gene_expression.csv",
                    )
                    if os.path.isfile(potential_dge_path):
                        dge_file_paths.append(potential_dge_path)
                    else:
                        print(
                            (
                                f"\t\t{potential_dge_path} does not exist. This means not enough cells "
                                "/ replicates were present to correctly assess differential gene expression."
                            ),
                            flush=True,
                        )
                current_cluster = node["cluster_id"]
                break
        current_noc = cluster_entry["number_of_clusters"]
        child_cluster_entry = cluster_entry

    # Orders the splits from root to leaf.
    dge_file_paths.reverse()
    cluster_size_weightings.reverse()
    return (dge_file_paths, cluster_size_weightings)


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
                    dge_files, cluster_size_weightings = get_branching_dges(
                        number_of_clusters=number_of_clusters,
                        cluster=cluster_id,
                        cluster_tree=relation_tree,
                        sub_folder=relative_sub_directory,
                    )

                    if len(dge_files) == 0:
                        print("\t\tNo DGE files have been detected. Skipping...", flush=True)
                        continue

                    print(
                        f"\t\tThe following DGE files have been selected: {dge_files}"
                    )
                    print(
                        f"\t\tThe following population size weights will be applied: {cluster_size_weightings}"
                    )
                    print("\t\tMerging gene expression data...")
                    merged_dge_file = load_and_merge_dges(
                        dge_paths=dge_files,
                        cluster_size_weights=cluster_size_weightings,
                    )
                    final_output_folder = os.path.join(
                        output_folder_path,
                        pathvalidate.sanitize_filename(str(number_of_clusters)),
                    )
                    print("\t\tWriting output to file...")
                    os.makedirs(final_output_folder, exist_ok=True)
                    final_output_file = os.path.join(
                        final_output_folder,
                        pathvalidate.sanitize_filename(
                            f"{cluster_id}_merged_dge_raw.csv"
                        ),
                    )
                    merged_dge_file.to_csv(
                        final_output_file,
                        sep=",",
                        encoding="utf-8",
                    )
                    print("\t\tWriting filtered and sorted output to file...")
                    final_output_file_filtered = os.path.join(
                        final_output_folder,
                        pathvalidate.sanitize_filename(
                            f"{cluster_id}_merged_dge_filtered.csv"
                        ),
                    )
                    merged_dge_file_filtered = merged_dge_file[
                        merged_dge_file[KEY_SIGNIFICANCE_DEPTH] >= 0
                    ].sort_values(
                        by=KEY_ENRICHMENT_SCORE_WEIGHTED, key=abs, ascending=False
                    )
                    merged_dge_file_filtered.to_csv(
                        final_output_file_filtered,
                        sep=",",
                        encoding="utf-8",
                    )
