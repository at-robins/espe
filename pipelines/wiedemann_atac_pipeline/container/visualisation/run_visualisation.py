#!/usr/bin/python
"""This module visualises peak data."""

import json
import os
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
MAX_RELEVANT_GENES_PER_DA = 100

print("Loading count matrix...", flush=True)
count_matrix_file = os.path.join(INPUT_FOLDER, "count_matrix_regularised.csv")
count_matrix = pd.read_csv(
    count_matrix_file, sep=",", header=0, index_col=0, encoding="utf-8"
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
plot = sns.clustermap(
    count_matrix.loc[sorted(relevant_peaks)],
    cmap="vlag",
    # z_score=0,
    standard_scale=0,
    center=0.5,
)
plot.savefig(
    os.path.join(
        MOUNT_PATHS["output"],
        "clustermap.svg",
    )
)
plt.close()
