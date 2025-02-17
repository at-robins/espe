#!/usr/bin/python
"""This module visualises motif-gene-associations."""

import csv
import json
import numpy as np
import os
import pandas as pd
import seaborn as sns

from pathlib import Path


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_MOTIF_DISCOVERY = MOUNT_PATHS["dependencies"]["motif_discovery"]
INPUT_FOLDER_ANNOTATIONS = MOUNT_PATHS["dependencies"]["peak_annotation"]
ANNOTATION_PATH = os.path.join(INPUT_FOLDER_ANNOTATIONS, "merged.mergedPeak")
MOTIF_FILE_NAME = "motif_annotated_peaks.tsv"
KEY_MOTIF_NAME = "Motif Name"
KEY_MOTIF_PEAK_ID = "PositionID"
KEY_ANNOTATION_GENE_NAME = "Gene Name"


def filter_genes_on_motif(group, filter_motif):
    """
    Returns true if the group contains the specified motif.
    """
    return filter_motif in group[KEY_MOTIF_NAME].values


def plot_correlation_clustermap(cm_data, cm_output_path, cluster_rows, cluster_columns):
    """
    Plots the correlation clustermaps.
    """
    # The diagonal might have extreme values, so median is better than mean.
    matrix_median = np.median(cm_data.values)
    # We want the actual standard deveation as we are looking at the whole population.
    matrix_sd = np.std(cm_data.values, ddof=0)
    upper_colour_cutoff = matrix_median + 3 * matrix_sd
    plot = sns.clustermap(
        cm_data,
        cbar_kws={"label": "motif correlation"},
        cmap="vlag",
        vmin=0,
        vmax=upper_colour_cutoff,
        center=1,
        row_cluster=cluster_rows,
        col_cluster=cluster_columns,
        figsize=(
            6.0 + 0.6 * len(cm_data.columns),
            4.0 + 0.3 * len(cm_data.index),
        ),
    )
    plot.ax_row_dendrogram.set_visible(False)
    plot.ax_col_dendrogram.set_visible(False)
    plot.ax_heatmap.grid(False)
    plot.savefig(cm_output_path)


print("Loading annotations...", flush=True)
annotation = pd.read_csv(
    ANNOTATION_PATH, sep="\t", header=0, index_col=0, encoding="utf-8"
)

for root, dirs, files in os.walk(INPUT_FOLDER_MOTIF_DISCOVERY):
    for file in files:
        if file == MOTIF_FILE_NAME:
            raw_path = Path(os.path.join(root, file))
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER_MOTIF_DISCOVERY)),
            )
            os.makedirs(output_folder_path, exist_ok=True)
            annotated_path = Path(
                os.path.join(output_folder_path, "annotated_motif_association.csv")
            )
            correlation_path_abs = Path(
                os.path.join(
                    output_folder_path, "motif_correlation_matrix_absolute.csv"
                )
            )
            correlation_path_rel = Path(
                os.path.join(
                    output_folder_path, "motif_correlation_matrix_relative.csv"
                )
            )
            clustermap_path_none = Path(
                os.path.join(
                    output_folder_path, "motif_correlation_relative_clustermap_none.svg"
                )
            )
            clustermap_path_both = Path(
                os.path.join(
                    output_folder_path, "motif_correlation_relative_clustermap_both.svg"
                )
            )
            clustermap_path_row = Path(
                os.path.join(
                    output_folder_path, "motif_correlation_relative_clustermap_row.svg"
                )
            )
            clustermap_path_col = Path(
                os.path.join(
                    output_folder_path,
                    "motif_correlation_relative_clustermap_column.svg",
                )
            )
            print(f"Processing {raw_path}...")
            print("\tAnnotating motifs...", flush=True)
            motif_table = pd.read_csv(
                raw_path, sep="\t", header=0, index_col=False, encoding="utf-8"
            )
            merged_motifs = pd.merge(
                left=motif_table,
                right=annotation[KEY_ANNOTATION_GENE_NAME],
                left_on=KEY_MOTIF_PEAK_ID,
                right_index=True,
                how="left",
            ).sort_values(KEY_ANNOTATION_GENE_NAME, ascending=True)
            print("\tSaving annotated motifs...", flush=True)
            merged_motifs.to_csv(
                annotated_path,
                sep=",",
                encoding="utf-8",
                index=False,
            )
            print("\tCreating motif correlation matrices...", flush=True)
            merged_motifs.dropna(
                axis=0, inplace=True, subset=[KEY_ANNOTATION_GENE_NAME, KEY_MOTIF_NAME]
            )
            genes = np.unique(merged_motifs[KEY_ANNOTATION_GENE_NAME].to_numpy())
            motifs = np.unique(merged_motifs[KEY_MOTIF_NAME].to_numpy())
            genes_grouped = merged_motifs.groupby(KEY_ANNOTATION_GENE_NAME)
            motif_occurences = merged_motifs.groupby(KEY_MOTIF_NAME).size()
            motif_background_frequency = motif_occurences / len(genes)
            # Creates zero filled data frames to write motif correlations to.
            motif_correlation_df_abs = pd.DataFrame(0, index=motifs, columns=motifs)
            motif_correlation_df_rel = pd.DataFrame(0, index=motifs, columns=motifs)
            for motif_a in motifs:
                genes_with_motif_a = genes_grouped.filter(
                    filter_genes_on_motif, True, *[motif_a]
                )
                motif_filtered_genes = len(
                    np.unique(genes_with_motif_a[KEY_ANNOTATION_GENE_NAME].to_numpy())
                )
                motif_filtered_motif_occurences = genes_with_motif_a.groupby(
                    KEY_MOTIF_NAME
                ).size()
                motif_correlation_df_abs[motif_a] = (
                    motif_filtered_motif_occurences / motif_filtered_genes
                )
                motif_correlation_df_rel[motif_a] = (
                    motif_correlation_df_abs[motif_a] / motif_background_frequency
                )

            # Empty values are the expected results of motifs not being present
            # in a motif subset, so these values should be filled with zeros.
            motif_correlation_df_abs.fillna(0, inplace=True)
            motif_correlation_df_rel.fillna(0, inplace=True)

            print("\tSaving motif correlation matrices...", flush=True)
            # Saves absolute correlation values.
            motif_correlation_df_abs.to_csv(
                correlation_path_abs,
                sep=",",
                encoding="utf-8",
                index=True,
            )
            # Saves correlation values normalised by per-gene motif background frequency.
            motif_correlation_df_rel.to_csv(
                correlation_path_rel,
                sep=",",
                encoding="utf-8",
                index=True,
            )

            print("\tPlotting motif correlation matrices...", flush=True)
            plot_correlation_clustermap(
                cm_data=motif_correlation_df_rel,
                cm_output_path=clustermap_path_none,
                cluster_rows=False,
                cluster_columns=False,
            )
            plot_correlation_clustermap(
                cm_data=motif_correlation_df_rel,
                cm_output_path=clustermap_path_row,
                cluster_rows=True,
                cluster_columns=False,
            )
            plot_correlation_clustermap(
                cm_data=motif_correlation_df_rel,
                cm_output_path=clustermap_path_col,
                cluster_rows=False,
                cluster_columns=True,
            )
            plot_correlation_clustermap(
                cm_data=motif_correlation_df_rel,
                cm_output_path=clustermap_path_both,
                cluster_rows=True,
                cluster_columns=True,
            )
