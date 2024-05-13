#!/usr/bin/python
"""This module plots marker genes."""

import anndata
import csv
import json
import os
import scanpy as sc

from matplotlib import pyplot as plt
from pathlib import PurePath

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
# CLUSTER_KEY = "number_of_clusters_17"


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


# print("Searching for data...")
# # Iterates over all sample directories and processes them conserving the directory structure.
# for root, dirs, files in os.walk(INPUT_FOLDER):
#     for file in files:
#         if file.casefold().endswith("cluster_relation.h5ad"):
#             file_path_adata = os.path.join(root, file)

#             output_folder_path = os.path.join(
#                 MOUNT_PATHS["output"],
#                 os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
#             )
#             os.makedirs(output_folder_path, exist_ok=True)
#             print(
#                 f"Processing file {file_path_adata}...",
#                 flush=True,
#             )
#             print(
#                 "\tReading data...",
#                 flush=True,
#             )
#             adata = anndata.read_h5ad(file_path_adata)
#             adata.obs[CLUSTER_KEY] = adata.obs[CLUSTER_KEY].astype("str")
#             print(
#                 "\tRanking genes...",
#                 flush=True,
#             )
#             sc.tl.rank_genes_groups(
#                 adata,
#                 groupby=CLUSTER_KEY,
#                 method="wilcoxon",
#             )
#             print(
#                 "\tPlotting...",
#                 flush=True,
#             )
#             fig = sc.pl.rank_genes_groups_stacked_violin(
#                 adata, n_genes=5, groupby=CLUSTER_KEY, return_fig=True
#             )
#             fig.savefig(
#                 os.path.join(output_folder_path, "marker_genes_violin.svg"),
#                 format="svg",
#             )


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


for directory_path in directory_paths:
    print(f"Processing directory {directory_path}...", flush=True)
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
        print(f"\t\tWriting total differentially expressed genes to file...", flush=True)
        with open(
            os.path.join(
                output_folder_path, "differentially_expressed_genes_total.json"
            ),
            "wt",
        ) as gene_file:
            json.dump(sorted(dge_set), gene_file)
