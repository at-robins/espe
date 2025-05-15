#!/usr/bin/python
"""This module performs gene set enrichment analysis."""

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
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
NETWORK_PATH_HUMAN_PROGENY = "/pathway_networks/human.pkl"
NETWORK_PATH_MOUSE_PROGENY = "/pathway_networks/mouse.pkl"
SIGNIFICANCE_THRESHOLD = 0.05
KEY_SCORE = "score"
KEY_SCORE_ABSOLUT = "absolute score"
KEY_PVALUE = "pvalue"
KEY_SIGNIFICANT = "significance"
VALUE_SIGNIFICANT_YES_UP = "more accessible"
VALUE_SIGNIFICANT_YES_DOWN = "less accessible"
VALUE_SIGNIFICANT_NO = "not significant"
KEY_PATHWAY_DB_GENESYMBOL = "target"


def parse_da_csv(dge_csv_path) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    # Only keeps relevant columns.
    annotated_dars = pd.read_csv(
        dge_csv_path, sep=",", header=0, index_col=None, encoding="utf-8"
    )[["padj", "log2FoldChange", "Gene Name"]].copy()
    annotated_dars.sort_values("padj", ascending=True, inplace=True)
    total_peaks = len(annotated_dars.index)

    # Only keeps relevant columns and renames the remaining ones.
    annotated_dars.drop(labels=["padj"], axis=1, inplace=True)
    annotated_dars.rename(columns={"log2FoldChange": KEY_SCORE}, inplace=True)

    # Only adds genes once for the most significant peak.
    annotated_dars.drop_duplicates(
        subset="Gene Name",
        keep="first",
        inplace=True,
    )

    # Remove invalid values.
    annotated_dars = annotated_dars[np.isfinite(annotated_dars[KEY_SCORE])].copy()

    # Set gene names as indices for the GSEA algorithm.
    annotated_dars.set_index("Gene Name", inplace=True)

    print(f"\tUsing {len(annotated_dars.index)} / {total_peaks} peaks.", flush=True)
    annotated_dars.sort_values(KEY_SCORE, ascending=False, inplace=True)

    return annotated_dars


def parse_da_csv_for_extension(dge_csv_path) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    # Only keeps relevant columns.
    annotated_dars = pd.read_csv(
        dge_csv_path, sep=",", header=0, index_col=None, encoding="utf-8"
    )[["log2FoldChange", "pvalue", "padj", "PeakID", "Gene Name"]].sort_values(
        "padj", ascending=True, inplace=False
    )

    # Only adds genes once for the most significant peak.
    annotated_dars.drop_duplicates(
        subset="Gene Name",
        keep="first",
        inplace=True,
    )

    # Set gene names as indices for the GSEA algorithm.
    annotated_dars.set_index("Gene Name", inplace=True)

    return annotated_dars


def row_to_significance(da_row):
    """
    Transforms a pvalue to a significance label.
    """
    if da_row[KEY_PVALUE] <= SIGNIFICANCE_THRESHOLD:
        if da_row[KEY_SCORE] >= 0:
            return VALUE_SIGNIFICANT_YES_UP
        else:
            return VALUE_SIGNIFICANT_YES_DOWN
    else:
        return VALUE_SIGNIFICANT_NO


def perform_enrichment_analysis(
    dge_path,
    pathway_db,
    output_folder,
    suffix,
):
    """
    Performs the pathway enrichment analysis.
    """
    os.makedirs(output_folder, exist_ok=True)

    print("\tParsing data data...")
    cluster_genes = parse_da_csv(dge_csv_path=dge_path)

    print("\tRunning GSEA...", flush=True)
    estimate, pvals = decoupler.run_mlm(
        mat=cluster_genes.T,
        net=pathway_db,
    )

    print("\tSaving results...", flush=True)
    pvals.rename(index={KEY_SCORE: KEY_PVALUE}, inplace=True)
    mlm_results = pd.concat([estimate, pvals]).T
    mlm_results.index.name = "pathway"
    mlm_results.sort_values(KEY_PVALUE, ascending=True, inplace=True)
    mlm_results.to_csv(
        os.path.join(
            output_folder,
            pathvalidate.sanitize_filename(f"pathway_analysis_{suffix}_table.csv"),
        ),
        sep=",",
        encoding="utf-8",
    )

    print("\tWriting per gene data...", flush=True)
    dge_lookup = parse_da_csv_for_extension(dge_csv_path=dge_path)
    pd.merge(
        pathway_db,
        dge_lookup,
        left_on=KEY_PATHWAY_DB_GENESYMBOL,
        right_index=True,
        how="left",
    ).to_csv(
        os.path.join(
            output_folder,
            pathvalidate.sanitize_filename(f"pathway_analysis_{suffix}_genes.csv"),
        ),
        sep=",",
        encoding="utf-8",
        index=False,
    )

    print("\tPlotting data...", flush=True)
    mlm_results.sort_values(KEY_SCORE, ascending=False, inplace=True)
    mlm_results[KEY_SIGNIFICANT] = pd.Categorical(
        mlm_results.apply(row_to_significance, axis=1),
        categories=[
            VALUE_SIGNIFICANT_YES_UP,
            VALUE_SIGNIFICANT_YES_DOWN,
            VALUE_SIGNIFICANT_NO,
        ],
        ordered=True,
    )
    mlm_results[KEY_SCORE_ABSOLUT] = np.abs(mlm_results[KEY_SCORE])

    height = 6
    width = 8
    fig, ax = plt.subplots(figsize=(width, height))
    colour_palette = sns.color_palette(["#E64B35", "#4DBBD5", "#B3B3B3FF"])
    sns.barplot(
        data=mlm_results,
        x=KEY_SCORE_ABSOLUT,
        y=mlm_results.index,
        hue=KEY_SIGNIFICANT,
        dodge=False,
        orient="h",
        width=0.75,
        palette=colour_palette,
        saturation=1.0,
        ax=ax,
    )
    # Adds a central line if positve and negative values are present.
    if np.min(mlm_results[KEY_SCORE]) < 0.0 and np.max(mlm_results[KEY_SCORE]) > 0.0:
        max_index_down = np.where(mlm_results[KEY_SCORE] < 0)[0][0]
        ax.axhline(y=max_index_down - 0.5, color="#000000FF", linestyle=":")
    ax.set(xlabel="absolute enrichment score", ylabel=None)
    legend = ax.get_legend()
    legend.set_title("Regulation")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            output_folder,
            pathvalidate.sanitize_filename(f"pathway_analysis_{suffix}_barplot.svg"),
        )
    )
    plt.close(fig)


print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    pathway_database_progeny = pd.read_pickle(NETWORK_PATH_HUMAN_PROGENY)
else:
    print("\tUsing mouse data...", flush=True)
    pathway_database_progeny = pd.read_pickle(NETWORK_PATH_MOUSE_PROGENY)

# Iterates over all directories and loads sample information.
print("Processing differential accessibility data...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.startswith("annotated_differential_accessibility") and file.endswith(
            ".csv"
        ):
            dge_file = os.path.join(root, file)
            print(f"Processing {dge_file}...", flush=True)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                file.removeprefix("annotated_differential_accessibility_").removesuffix(
                    ".csv"
                ),
            )
            print("\tUsing progeny database...", flush=True)
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_progeny,
                output_folder=output_folder_path,
                suffix="progeny",
            )
