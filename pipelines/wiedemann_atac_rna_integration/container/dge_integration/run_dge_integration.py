#!/usr/bin/python
"""This module integrates ATAC and RNA-Seq data."""

import csv
import json
import numpy as np
import os
import pandas as pd
import seaborn as sns

from collections import OrderedDict
from matplotlib.patches import Patch
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "monospace"

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["input"]
KEY_GENE_SYMBOL = "Gene symbol"
KEY_LFC = "log2FC"
KEY_SUFFIX_RNA = " RNA"
KEY_SUFFIX_ATAC = " ATAC"
ENRICHMENT_CATEGORY_CUTOFF = 1.0
KEY_ENRICHMENT = "enrichment"
ENRICHMENT_VALUE_TRUE = "log₂ fold changes ≥ 1.0"
ENRICHMENT_VALUE_FALSE = "log₂ fold changes < 1.0"


def parse_csv_file(csv_path) -> pd.DataFrame:
    """
    Parses the DGE / DA CSV file.
    """
    csv_table = pd.read_csv(
        csv_path,
        sep=",",
        header=0,
        index_col=None,
        encoding="utf-8",
    )
    return csv_table[csv_table["FDR"] <= 0.05].copy()


def lfcs_to_category(da_row):
    """
    Converts log2-fold-changes to plotable enrichment categories.
    """
    lfc_rna = abs(float(da_row[KEY_LFC + KEY_SUFFIX_RNA]))
    lfc_atac = abs(float(da_row[KEY_LFC + KEY_SUFFIX_ATAC]))
    return (
        ENRICHMENT_VALUE_TRUE
        if lfc_atac >= ENRICHMENT_CATEGORY_CUTOFF
        and lfc_rna >= ENRICHMENT_CATEGORY_CUTOFF
        else ENRICHMENT_VALUE_FALSE
    )


print(
    "Processing differential accessibility and differential gene expression data...",
    flush=True,
)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file == "atac_da.csv":
            da_file = os.path.join(root, file)
            dge_file = os.path.join(root, "rna_dge.csv")
            print(f"Processing {da_file} and {dge_file}...", flush=True)
            directory_path_output = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            os.makedirs(directory_path_output, exist_ok=True)

            # Loads and prefilters the data.
            da_table = parse_csv_file(da_file)
            dge_table = parse_csv_file(dge_file)

            print("\tMerging data...", flush=True)
            merged_table = pd.merge(
                da_table,
                dge_table,
                on=KEY_GENE_SYMBOL,
                how="inner",
                suffixes=[KEY_SUFFIX_ATAC, KEY_SUFFIX_RNA],
            )

            reduced_accessibility_table = merged_table[
                merged_table[KEY_LFC + KEY_SUFFIX_ATAC] < 0.0
            ]
            increased_accessibility_table = merged_table[
                merged_table[KEY_LFC + KEY_SUFFIX_ATAC] > 0.0
            ]
            # Deduplicates peaks assigned to the same gene
            # by using maximum / minimum values.
            reduced_accessibility_indices_deduplicated = (
                reduced_accessibility_table.groupby(KEY_GENE_SYMBOL)[
                    KEY_LFC + KEY_SUFFIX_ATAC
                ].idxmin()
            )

            increased_accessibility_indices_deduplicated = (
                increased_accessibility_table.groupby(KEY_GENE_SYMBOL)[
                    KEY_LFC + KEY_SUFFIX_ATAC
                ].idxmax()
            )
            merged_deduplicated_indices = []
            merged_deduplicated_indices.extend(
                reduced_accessibility_indices_deduplicated
            )
            merged_deduplicated_indices.extend(
                increased_accessibility_indices_deduplicated
            )

            merged_table = merged_table.loc[merged_deduplicated_indices].copy()
            # Skips samples if no overlapping features.
            overlapping_features = len(merged_table.index)
            if overlapping_features > 0:
                print(
                    f"\tDetected {overlapping_features} overlapping features...",
                    flush=True,
                )
            else:
                print("\tNo overlap detected. Skipping samples...", flush=True)
                continue

            # Sets data point colour variable.
            merged_table[KEY_ENRICHMENT] = pd.Categorical(
                merged_table.apply(lfcs_to_category, axis=1),
                categories=[ENRICHMENT_VALUE_TRUE, ENRICHMENT_VALUE_FALSE],
                ordered=True,
            )

            print("\tPlotting data...", flush=True)
            lfc_axis_scale_rna = (
                merged_table[KEY_LFC + KEY_SUFFIX_RNA].abs().max() * 1.1
            )
            lfc_axis_scale_atac = (
                merged_table[KEY_LFC + KEY_SUFFIX_ATAC].abs().max() * 1.1
            )
            fig, ax = plt.subplots(figsize=(6, 6))
            sns.scatterplot(
                data=merged_table,
                x=KEY_LFC + KEY_SUFFIX_ATAC,
                y=KEY_LFC + KEY_SUFFIX_RNA,
                hue=KEY_ENRICHMENT,
                palette={
                    ENRICHMENT_VALUE_TRUE: "#000000FF",
                    ENRICHMENT_VALUE_FALSE: "#B3B3B3FF",
                },
                ax=ax,
            )
            # Set axes limits.
            ax.set_xlim(-lfc_axis_scale_atac, lfc_axis_scale_atac)
            ax.set_ylim(-lfc_axis_scale_rna, lfc_axis_scale_rna)
            # Adds quadrant lines.
            ax.axvline(x=0, color="#000000FF", linestyle=":")
            ax.axhline(y=0, color="#000000FF", linestyle=":")
            # Sets labels.
            ax.set(
                xlabel="genomic accessibility log₂ fold change",
                ylabel="transcript abundance log₂ fold change",
            )
            legend = ax.get_legend()
            legend.set_title("Enrichment")

            # Adds labels for enriched features.
            def label_enriched_features(df_row, plot_axis):
                """
                Local label function with access to the axis variable.
                """
                if df_row[KEY_ENRICHMENT] == ENRICHMENT_VALUE_TRUE:
                    plot_axis.text(
                        df_row[KEY_LFC + KEY_SUFFIX_ATAC] + 0.02 * lfc_axis_scale_atac,
                        df_row[KEY_LFC + KEY_SUFFIX_RNA] + 0.005 * lfc_axis_scale_rna,
                        df_row[KEY_GENE_SYMBOL],
                    )

            merged_table.apply(label_enriched_features, args=(ax,), axis=1)
            # Saves and closes the plot.
            fig.tight_layout()
            fig.savefig(
                os.path.join(
                    directory_path_output, "atac_rna_integration_scatterplot.svg"
                )
            )
            plt.close(fig)
