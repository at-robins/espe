#!/usr/bin/python
"""This module visualises motif-gene-associations."""

import csv
import json
import numpy as np
import os
import pandas as pd

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
    return filter_motif in group[KEY_MOTIF_NAME].values


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
                os.path.join(output_folder_path, f"annotated_motif_association.csv")
            )
            correlation_path_abs = Path(
                os.path.join(
                    output_folder_path, f"motif_correlation_matrix_absolute.csv"
                )
            )
            correlation_path_rel = Path(
                os.path.join(
                    output_folder_path, f"motif_correlation_matrix_relative.csv"
                )
            )
            print(f"Annotating {raw_path}...", flush=True)
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
            # Writes annotated motifs to CSV.
            merged_motifs.to_csv(
                annotated_path,
                sep=",",
                encoding="utf-8",
                index=False,
            )
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
