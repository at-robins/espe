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
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["cluster_relation_dge"] + "/"
GMT_PATH_HUMAN_REACTOME = "/gmts/reactome_human.gmt"
GMT_PATH_MOUSE_REACTOME = "/gmts/reactome_mouse.gmt"
GMT_PATH_HUMAN_GO = "/gmts/go_human.gmt"
GMT_PATH_MOUSE_GO = "/gmts/go_mouse.gmt"
KEY_PATHWAY_DB_GENESET = "geneset"
KEY_PATHWAY_DB_GENESYMBOL = "genesymbol"
FILTER_PATHWAY_DB_MIN = 15
FILTER_PATHWAY_DB_MAX = 500
PLOT_PATHWAYS_MAX = 30


def parse_gmt(pth: Path) -> pd.DataFrame:
    """
    Parses a GMT file as decoupler pathway dataframe.
    """
    pathways = {}

    with Path(pth).open("r", encoding="utf-8") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    df = pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=[KEY_PATHWAY_DB_GENESET, KEY_PATHWAY_DB_GENESYMBOL],
    )

    # Filters the pathway database by size.
    geneset_size = df.groupby(KEY_PATHWAY_DB_GENESET).size()
    gsea_genesets = geneset_size.index[
        (geneset_size > FILTER_PATHWAY_DB_MIN) & (geneset_size < FILTER_PATHWAY_DB_MAX)
    ]
    df = df[df[KEY_PATHWAY_DB_GENESET].isin(gsea_genesets)]

    return df


def parse_dge_csv(dge_csv_path) -> pd.DataFrame:
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
            f_scores.append(float(row["F"]) * np.sign(float(row["logFC"])))
        df = pd.DataFrame(data={"scores": f_scores}, index=genes)
        df.sort_values(by=["scores"], inplace=True, ascending=False)

        return df


def parse_dge_csv_for_extension(dge_csv_path) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    with open(dge_csv_path, newline="", encoding="utf-8") as csvfile:
        dge_reader = csv.DictReader(
            csvfile, dialect="unix", delimiter=",", quotechar='"'
        )
        genes = []
        remainder = {
            "logFC": [],
            "PValue": [],
            "FDR": [],
        }
        for row in dge_reader:
            genes.append(row.pop(""))
            remainder["logFC"].append(row["logFC"])
            remainder["PValue"].append(row["PValue"])
            remainder["FDR"].append(row["FDR"])
        df = pd.DataFrame(data=remainder, index=genes)

        return df


def perform_enrichment_analysis(dge_path, pathway_db, output_folder, suffix):
    """
    Performs the pathway enrichment analysis.
    """
    os.makedirs(output_folder, exist_ok=True)

    cluster_genes = parse_dge_csv(
        dge_csv_path=dge_path,
    )
    scores, norm, pvals = decoupler.run_gsea(
        cluster_genes.T,
        pathway_db,
        source=KEY_PATHWAY_DB_GENESET,
        target=KEY_PATHWAY_DB_GENESYMBOL,
    )

    gsea_results = (
        pd.concat({"score": scores.T, "norm": norm.T, "pval": pvals.T}, axis=1)
        .droplevel(level=1, axis=1)
        .assign(
            **{"abs_norm": lambda x: np.abs(x["norm"])},
            **{
                "regulation": lambda x: pd.Categorical(
                    map(
                        lambda z: ("down-regulated" if z < 0.0 else "up-regulated"),
                        np.sign(x["score"]),
                    ),
                    categories=["up-regulated", "down-regulated"],
                    ordered=True,
                )
            },
        )
        .sort_values(["pval", "abs_norm"])
    )
    gsea_results.to_csv(
        os.path.join(
            output_folder, pathvalidate.sanitize_filename(f"gsea_{suffix}_table.csv")
        ),
        sep=",",
        encoding="utf-8",
    )

    print("\tWriting per gene data...")
    gsea_results_filtered_sets = gsea_results[gsea_results["pval"] <= 0.05].sort_values(
        "norm", ascending=False
    ).index.values
    dge_lookup = parse_dge_csv_for_extension(
        dge_csv_path=dge_path,
    )
    pd.merge(
        pathway_db[pathway_db[KEY_PATHWAY_DB_GENESET].isin(gsea_results_filtered_sets)],
        dge_lookup,
        left_on=KEY_PATHWAY_DB_GENESYMBOL,
        right_index=True,
        how="left",
    ).to_csv(
        os.path.join(
            output_folder, pathvalidate.sanitize_filename(f"gsea_{suffix}_genes.csv")
        ),
        sep=",",
        encoding="utf-8",
    )

    print("\tPlotting data...")
    gsea_results_filtered = (
        gsea_results[gsea_results["pval"] <= 0.05]
        .head(30)
        .sort_values("norm", ascending=False)
    )

    if len(gsea_results_filtered) == 0:
        print("\tNo significantly altered pathways. Skipping sample...")
    else:
        height = (7 * (len(gsea_results_filtered) / PLOT_PATHWAYS_MAX)) + 1
        width = (max(map(len, gsea_results_filtered.index.values)) / 10) + 4
        fig, ax = plt.subplots(figsize=(15, height))
        colour_palette = sns.color_palette(["#E64B35", "#4DBBD5"])
        sns.barplot(
            data=gsea_results_filtered,
            x="abs_norm",
            y=gsea_results_filtered.index,
            hue="regulation",
            dodge=False,
            orient="h",
            width=0.75,
            palette=colour_palette,
            saturation=1.0,
            ax=ax,
        )
        ax.set(xlabel="absolute NES", ylabel=None)
        legend = ax.get_legend()
        legend.set_title("Regulation")
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        fig.tight_layout()
        fig.savefig(
            os.path.join(
                output_folder,
                pathvalidate.sanitize_filename(f"gsea_{suffix}_barplot.svg"),
            )
        )
        plt.close(fig)


print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    pathway_database_reactome = parse_gmt(GMT_PATH_HUMAN_REACTOME)
    pathway_database_go = parse_gmt(GMT_PATH_HUMAN_GO)
else:
    print("\tUsing mouse data...", flush=True)
    pathway_database_reactome = parse_gmt(GMT_PATH_MOUSE_REACTOME)
    pathway_database_go = parse_gmt(GMT_PATH_MOUSE_GO)

# Iterates over all directories and loads sample information.
print("Processing differential gene expression data...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith("differential_gene_expression.csv"):
            dge_file = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"], root.removeprefix(INPUT_FOLDER)
            )
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_go,
                output_folder=output_folder_path,
                suffix="go",
            )
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_reactome,
                output_folder=output_folder_path,
                suffix="rectome",
            )
