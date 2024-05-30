#!/usr/bin/python
"""This module visualises peak data."""

import csv
import json
import os
import pandas as pd
import seaborn as sns

from collections import OrderedDict
from matplotlib.patches import Patch
from matplotlib import pyplot as plt

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
MAX_RELEVANT_GENES_PER_DA = 100
COLOUR_PALETTE = [
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148",
    "#B09C85",
]


def count_matrix_sample_headers(count_matrix_file_path) -> [str]:
    """
    Returns the headers of the samples present in the count matrix.
    """
    with open(
        count_matrix_file_path, newline="", mode="rt", encoding="utf-8"
    ) as count_matrix_in:
        count_matrix_reader = csv.reader(
            count_matrix_in, dialect="unix", delimiter=",", quotechar='"'
        )
        matrix_row = next(count_matrix_reader)
        return matrix_row[1 : len(matrix_row)]


print("Loading count matrix...", flush=True)
count_matrix_file = os.path.join(INPUT_FOLDER, "count_matrix_regularised.csv")
count_matrix_conditions = count_matrix_sample_headers(count_matrix_file)

# Initialises the column name map with replicate 1 for each condition.
count_matrix_header_map = {}
for name in set(count_matrix_conditions):
    count_matrix_header_map[name] = 1

# Creates a unique column name per replicate.
count_matrix_col_names = []
for name in count_matrix_conditions:
    count_matrix_col_names.append(f"{name} ({count_matrix_header_map[name]})")
    count_matrix_header_map[name] = count_matrix_header_map[name] + 1

# Loads the actual count matrix.
count_matrix = pd.read_csv(
    count_matrix_file,
    sep=",",
    header=0,
    names=count_matrix_col_names,
    index_col=0,
    encoding="utf-8",
)

print("Loading relevant peaks...", flush=True)
relevant_peaks = set()
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.startswith("differential_accessibility") and file.endswith(".csv"):
            da_table = pd.read_csv(
                os.path.join(root, file),
                sep=",",
                header=0,
                index_col=0,
                encoding="utf-8",
            )
            relevant_peaks = relevant_peaks | set(
                da_table[da_table["padj"] >= 0.05]
                .sort_values("padj", ascending=True)
                .head(MAX_RELEVANT_GENES_PER_DA)
                .index
            )
print(f"Selected {len(relevant_peaks)} relevant peaks...", flush=True)

print("Plotting clustermap...", flush=True)
# The set is sorted into a list afterwards to allow consistent assignment
# colour assignment per label.
unique_col_names = sorted(set(count_matrix_conditions), key=str.casefold)
col_colour_palette = sns.color_palette(
    palette=COLOUR_PALETTE, n_colors=len(unique_col_names)
)
lut = OrderedDict(zip(unique_col_names, col_colour_palette))
# Assigns the column names colour based on their respective condition.
column_colours = {}
for index_col, name_col in enumerate(count_matrix_col_names):
    column_colours[name_col] = lut[count_matrix_conditions[index_col]]

plot = sns.clustermap(
    count_matrix.loc[sorted(relevant_peaks)],
    cmap="vlag",
    col_colors=pd.Series(column_colours, name="Sample type"),
    cbar_kws={"label": "standardised rlog counts"},
    standard_scale=0,
    center=0.5,
)
plot.ax_row_dendrogram.set_visible(False)

# Plots a legend for the per sample type / condition colour code.
handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(
    handles,
    lut,
    title="Sample type",
    bbox_to_anchor=(0.15, 0.5),
    bbox_transform=plt.gcf().transFigure,
    # Where the anchor is on the legend box.
    loc="center right",
)

plot.savefig(
    os.path.join(
        MOUNT_PATHS["output"],
        "clustermap.svg",
    )
)
plt.close()
