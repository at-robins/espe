#!/usr/bin/python
"""This module performs gene set enrichment analysis."""

import anndata
import csv

print("Loading decoupler...", flush=True)
import decoupler
import json
import numpy as np
import os
import pandas as pd
import pathvalidate
import seaborn as sns

from itertools import chain, repeat
from matplotlib import pyplot as plt
from pathlib import Path


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
FOLDER_COMPARISON_CLUSTER = "cluster_comparison"
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["cluster_relation_dge"] + "/"

GMT_PATH_HUMAN = "/gmts/reactome_human.gmt"
GMT_PATH_MOUSE = "/gmts/reactome_mouse.gmt"
KEY_PATHWAY_DB_GENESET = "geneset"
KEY_PATHWAY_DB_GENESYMBOL = "genesymbol"
FILTER_PATHWAY_DB_MIN = 15
FILTER_PATHWAY_DB_MAX = 500
KEY_CLUSTERING_GROUP = "clustering_group"
KEY_CLUSTER = "cluster"
KEY_DGE_FILE = "dge_file"


def parse_gmt(pth: Path) -> pd.DataFrame:
    """
    Parses a GMT file as decoupler pathway dataframe.
    """
    pathways = {}

    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=[KEY_PATHWAY_DB_GENESET, KEY_PATHWAY_DB_GENESYMBOL],
    )


def parse_dge_csv(dge_csv_path, output_folder) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        genes = []
        f_scores = []
        for row in dge_reader:
            genes.append(row[""])
            f_scores.append(float(row["F"]))
        df = pd.DataFrame(data={"scores": f_scores}, index=genes)
        df.sort_values(by=["scores"], inplace=True, ascending=False)
        return df


print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    pathway_database = parse_gmt(GMT_PATH_HUMAN)
else:
    print("\tUsing mouse data...", flush=True)
    pathway_database = parse_gmt(GMT_PATH_MOUSE)

# Filters the pathway database by size.
geneset_size = pathway_database.groupby(KEY_PATHWAY_DB_GENESET).size()
gsea_genesets = geneset_size.index[
    (geneset_size > FILTER_PATHWAY_DB_MIN) & (geneset_size < FILTER_PATHWAY_DB_MAX)
]
pathway_database = pathway_database[
    pathway_database[KEY_PATHWAY_DB_GENESET].isin(gsea_genesets)
]

# Iterates over all directories and loads sample information.
print("Loading cluster information...", flush=True)
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
    for cluster_info in group_info:
        print(f"\tProcessing cluster {cluster_info[KEY_CLUSTER]}...", flush=True)
        output_path = os.path.join(
            MOUNT_PATHS["output"],
            pathvalidate.sanitize_filename(clustering_group),
            pathvalidate.sanitize_filename(cluster_info[KEY_CLUSTER]),
        )
        os.makedirs(output_path, exist_ok=True)
        cluster_genes = parse_dge_csv(
            cluster_info[KEY_DGE_FILE],
            output_folder=output_path,
        )
        scores, norm, pvals = decoupler.run_gsea(
            cluster_genes.T,
            pathway_database,
            source=KEY_PATHWAY_DB_GENESET,
            target=KEY_PATHWAY_DB_GENESYMBOL,
        )

        gsea_results = (
            pd.concat({"score": scores.T, "norm": norm.T, "pval": pvals.T}, axis=1)
            .droplevel(level=1, axis=1)
            .sort_values("pval")
        )
        gsea_results.to_csv(
            os.path.join(output_path, "gsea_table.csv"), sep=",", encoding="utf-8"
        )

        print("\tPlotting data...")
        gsea_results_filtered = (
            gsea_results[gsea_results["pval"] <= 0.05]
            .head(20)
            .assign(
                **{"-log10(pval)": lambda x: -np.log10(x["pval"])},
                **{
                    "regulation": lambda x: pd.Categorical(
                        map(
                            lambda z: "down-regulated" if z < 0.0 else "up-regulated",
                            np.sign(x["score"]),
                        ),
                        categories=["down-regulated", "up-regulated"],
                        ordered=True,
                    )
                },
            )
        )
        fig, ax = plt.subplots(figsize=(15, 5))
        colour_palette = sns.color_palette(["#7497F5", "#EA7B60"])
        sns.barplot(
            data=gsea_results_filtered,
            x="-log10(pval)",
            y=gsea_results_filtered.index,
            hue="regulation",
            orient="h",
            palette=colour_palette,
            saturation=1.0,
            ax=ax,
        )
        ax.set(xlabel="-log₁₀(p-value)", ylabel=None)
        legend = ax.get_legend()
        legend.set_title("Regulation")
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        fig.tight_layout()
        fig.savefig(os.path.join(output_path, "gsea_barplot_pval.svg"))
        plt.close(fig)
