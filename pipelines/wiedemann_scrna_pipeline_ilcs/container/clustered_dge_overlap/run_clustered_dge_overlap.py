#!/usr/bin/python
"""This module calculates differentially expressed genes."""

import anndata
import csv
import json
import os
import pathvalidate


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
FOLDER_COMPARISON_CLUSTER = "cluster_comparison"
INPUT_FOLDER = (
    os.path.join(
        MOUNT_PATHS["dependencies"]["clustered_differential_gene_expression"],
        FOLDER_COMPARISON_CLUSTER,
    )
    + "/"
)
KEY_CLUSTERING_GROUP = "clustering_group"
KEY_CLUSTER = "cluster"
KEY_DGE_FILE = "dge_file"
FDR_CUTOFF = 0.05
LFC_CUTOFF = 0.5


def parse_dge_csv(dge_csv_path, cluster, gene_array):
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        for row in dge_reader:
            gene = row[""]
            fdr = float(row["FDR"])
            lfc = float(row["logFC"])
            if fdr <= FDR_CUTOFF and abs(lfc) >= LFC_CUTOFF:
                if lfc > 0:
                    regulation = "up-regulated"
                else:
                    regulation = "down-regulated"
                is_new_gene_regulation = True
                for gene_regulation_info in gene_array:
                    # Adds the current cluster if the same gene regulation has been identified in another cluster.
                    if (
                        gene_regulation_info["gene"] == gene
                        and gene_regulation_info["regulation"] == regulation
                    ):
                        gene_regulation_info["clusters"].append(cluster)
                        is_new_gene_regulation = False
                        break

                if is_new_gene_regulation:
                    gene_array.append(
                        {"gene": gene, "regulation": regulation, "clusters": [cluster]}
                    )


print("Loading cluster information...", flush=True)

# Iterates over all directories and loads sample information.
dge_input = {}
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("info.csv"):
            with open(
                os.path.join(root, file), newline="", encoding="utf-8"
            ) as csvfile:
                info_reader = csv.DictReader(
                    csvfile, dialect="unix", delimiter=",", quotechar='"'
                )
                first_row = next(info_reader)
                clustering_group = first_row[KEY_CLUSTERING_GROUP]
                cluster = first_row[KEY_CLUSTER]
                cluster_info = {
                    KEY_CLUSTER: cluster,
                    KEY_DGE_FILE: os.path.join(
                        root, "differential_gene_expression.csv"
                    ),
                }

                group_info = dge_input.get(clustering_group)
                if group_info is None:
                    dge_input[clustering_group] = [cluster_info]
                else:
                    group_info.append(cluster_info)

print("Processing differential gene expression data...", flush=True)
for clustering_group, group_info in dge_input.items():
    print(f"Processing clustering group {clustering_group}...", flush=True)
    gene_overlap = []
    for cluster_info in group_info:
        print(f"\tProcessing cluster {cluster_info[KEY_CLUSTER]}...", flush=True)
        parse_dge_csv(
            cluster_info[KEY_DGE_FILE],
            cluster=cluster_info[KEY_CLUSTER],
            gene_array=gene_overlap,
        )
    gene_overlap = sorted(
        gene_overlap,
        key=lambda gene_regulation_info: len(gene_regulation_info["clusters"]),
    )

    output_path = os.path.join(
        MOUNT_PATHS["output"],
        pathvalidate.sanitize_filename(clustering_group),
    )
    os.makedirs(output_path, exist_ok=True)
    with open(
        os.path.join(output_path, "overlapping_dge.csv"),
        mode="w",
        newline="",
        encoding="utf-8",
    ) as csvfile:
        info_writer = csv.writer(
            csvfile,
            dialect="unix",
            delimiter=",",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
        )
        info_writer.writerows(
            [
                [
                    "Gene",
                    "Regultaion",
                    "Overlap",
                    "Clusters",
                ],
            ]
        )
        for gene_regulation_info in gene_overlap:
            info_writer.writerows(
                [
                    [
                        gene_regulation_info["gene"],
                        gene_regulation_info["regulation"],
                        len(gene_regulation_info["clusters"]),
                        ";".join(sorted(gene_regulation_info["clusters"])),
                    ],
                ]
            )
