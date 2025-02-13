#!/usr/bin/python
"""This module visualises motif-gene-associations."""

import csv
import json
import os
import pandas as pd

from pathlib import Path


MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER_MOTIF_DISCOVERY = MOUNT_PATHS["dependencies"]["motif_discovery"]
INPUT_FOLDER_ANNOTATIONS = MOUNT_PATHS["dependencies"]["peak_annotation"]
ANNOTATION_PATH = os.path.join(INPUT_FOLDER_ANNOTATIONS, "merged.mergedPeak")
MOTIF_FILE_NAME = "motif_annotated_peaks.tsv"
KEY_MOTIF_PEAK_ID = "PositionID"
KEY_ANNOTATION_GENE_NAME = "Gene Name"


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
            annotated_path = Path(os.path.join(output_folder_path, f"annotated_{file}"))
            print(f"Annotating {raw_path}...", flush=True)
            motif_table = pd.read_csv(
                raw_path, sep="\t", header=0, index_col=False, encoding="utf-8"
            )
            pd.merge(
                left=motif_table,
                right=annotation[KEY_ANNOTATION_GENE_NAME],
                left_on=KEY_MOTIF_PEAK_ID,
                right_index=True,
                how="left",
            ).sort_values(KEY_ANNOTATION_GENE_NAME, ascending=True).to_csv(
                annotated_path,
                sep=",",
                encoding="utf-8",
                index=False,
            )
