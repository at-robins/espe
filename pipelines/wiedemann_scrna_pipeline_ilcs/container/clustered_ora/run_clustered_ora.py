#!/usr/bin/python
"""This module performs gene overrepresentation analysis."""

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
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["clustered_dge_overlap"] + "/"

GMT_PATH_HUMAN = "/gmts/reactome_human.gmt"
GMT_PATH_MOUSE = "/gmts/reactome_mouse.gmt"
KEY_PATHWAY_DB_GENESET = "geneset"
KEY_PATHWAY_DB_GENESYMBOL = "genesymbol"
FILTER_PATHWAY_DB_MIN = 15
FILTER_PATHWAY_DB_MAX = 500
FILTER_SAMPLE_MIN_GENES = 5
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


def parse_overlap_csv(overlap_csv_path) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(overlap_csv_path, newline="", encoding="utf-8") as csvfile:
        overlap_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        overlaps = {}
        for row in overlap_reader:
            # Only use genes unique for a cluster
            if int(row["Overlap"]) == 1:
                gene = row["Gene"]
                cluster = row["Clusters"]
                regulation = row["Regulation"]
                if overlaps.get(cluster) is None:
                    overlaps[cluster] = {regulation: [gene]}
                else:
                    if overlaps[cluster].get(regulation) is None:
                        overlaps[cluster][regulation] = [gene]
                    else:
                        overlaps[cluster][regulation].append(gene)

        return overlaps


def perform_ora(gene_set, pathway_db, output_folder, specifier):
    """
    Runs overrepresentation analysis.
    """
    if len(gene_set) < FILTER_SAMPLE_MIN_GENES:
        print(
            f"\tThe gene set {gene_set} does not contain enough genes for testing. Skipping sample...",
            flush=True,
        )
        return

    os.makedirs(output_folder, exist_ok=True)
    # gene_set_df = pd.DataFrame(data={"gene": np.repeat(1.0, len(gene_set))}, index=gene_set)
    ora_results = decoupler.get_ora_df(
        gene_set,
        pathway_database,
        source=KEY_PATHWAY_DB_GENESET,
        target=KEY_PATHWAY_DB_GENESYMBOL,
    )
    ora_results.sort_values("FDR p-value", inplace=True)

    ora_results.to_csv(
        os.path.join(output_folder, f"ora_table_{pathvalidate.sanitize_filename(specifier)}.csv"), sep=",", encoding="utf-8"
    )

    print("\tPlotting data...")
    ora_results_filtered = (
        ora_results[ora_results["FDR p-value"] <= 0.05]
        .head(20)
        .assign(**{"-log10(pval)": lambda x: -np.log10(x["FDR p-value"])})
    )
    if len(ora_results_filtered) < 1:
        print("\tNot enough data for plotting...")
        return

    fig, ax = plt.subplots(figsize=(15, 5))
    # colour_palette = sns.color_palette(["#7497F5", "#EA7B60"])
    sns.barplot(
        data=ora_results_filtered,
        x="-log10(pval)",
        y="Term",
        # hue=None,
        orient="h",
        # palette=colour_palette,
        saturation=1.0,
        ax=ax,
    )
    ax.set(xlabel="-log₁₀(FDR)", ylabel=None)
    # legend = ax.get_legend()
    # legend.set_title("Regulation")
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            output_folder,
            f"ora_barplot_fdr_{pathvalidate.sanitize_filename(specifier)}.svg",
        )
    )
    plt.close(fig)


print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and skip_env == "human":
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
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("overlapping_dge.csv"):
            output_base_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            overlap_input_path = os.path.join(root, file)
            print(
                f"Running overrepresentation analysis for sample {overlap_input_path}...",
                flush=True,
            )
            for cluster, regulation_data in parse_overlap_csv(
                overlap_input_path
            ).items():
                print(f"\tRunning overrepresentation analysis for cluster {cluster}...")
                all_genes = []
                for regulation, genes in regulation_data.items():
                    all_genes.extend(genes.copy())
                    perform_ora(
                        gene_set=genes,
                        pathway_db=pathway_database,
                        output_folder=os.path.join(
                            output_base_path, pathvalidate.sanitize_filename(cluster)
                        ),
                        specifier=regulation,
                    )
                perform_ora(
                    gene_set=all_genes,
                    pathway_db=pathway_database,
                    output_folder=os.path.join(
                        output_base_path, pathvalidate.sanitize_filename(cluster)
                    ),
                    specifier="all",
                )
