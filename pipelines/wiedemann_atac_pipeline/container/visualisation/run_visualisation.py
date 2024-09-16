#!/usr/bin/python
"""This module visualises peak data."""

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

ANNOTATION_RENAMING_MAP = OrderedDict(
    [
        ("exon", "exon"),
        ("intron", "intron"),
        ("Intergenic", "intergenic"),
        ("non-coding", "non-coding"),
        ("promoter-TSS", "promoter"),
        ("5' UTR", "5' UTR"),
        ("3' UTR", "3' UTR"),
        ("TTS", "TTS"),
    ]
)

LOCATION_PLOT_LABEL_POS_HIGH = "       LFC ≥  1.0"
LOCATION_PLOT_LABEL_POS_MID = " 1.0 > LFC ≥  0.5"
LOCATION_PLOT_LABEL_POS_LOW = " 0.5 > LFC ≥  0.0"
LOCATION_PLOT_LABEL_NEG_LOW = " 0.0 > LFC ≥ -0.5"
LOCATION_PLOT_LABEL_NEG_MID = "-0.5 > LFC ≥ -1.0"
LOCATION_PLOT_LABEL_NEG_HIGH = "-1.0 > LFC       "

LOCATION_PLOT_COLOURS = {
    LOCATION_PLOT_LABEL_POS_HIGH: "#E64B35",
    LOCATION_PLOT_LABEL_POS_MID: "#EE8778",
    LOCATION_PLOT_LABEL_POS_LOW: "#F7C3BC",
    LOCATION_PLOT_LABEL_NEG_LOW: "#C4E8F1",
    LOCATION_PLOT_LABEL_NEG_MID: "#88D2E3",
    LOCATION_PLOT_LABEL_NEG_HIGH: "#4DBBD5",
}

LOCATION_PLOT_BINS = OrderedDict(
    [
        (LOCATION_PLOT_LABEL_POS_HIGH, 1.0),
        (LOCATION_PLOT_LABEL_POS_MID, 0.5),
        (LOCATION_PLOT_LABEL_POS_LOW, 0.0),
        (LOCATION_PLOT_LABEL_NEG_LOW, -0.5),
        (LOCATION_PLOT_LABEL_NEG_MID, -1.0),
        (LOCATION_PLOT_LABEL_NEG_HIGH, float("-inf")),
    ]
)

LOCATION_PLOT_LABEL_ORDER_STACK = [
    LOCATION_PLOT_LABEL_POS_LOW,
    LOCATION_PLOT_LABEL_POS_MID,
    LOCATION_PLOT_LABEL_POS_HIGH,
    LOCATION_PLOT_LABEL_NEG_LOW,
    LOCATION_PLOT_LABEL_NEG_MID,
    LOCATION_PLOT_LABEL_NEG_HIGH,
]
LOCATION_PLOT_LABEL_ORDER_LEGEND = [
    LOCATION_PLOT_LABEL_POS_HIGH,
    LOCATION_PLOT_LABEL_POS_MID,
    LOCATION_PLOT_LABEL_POS_LOW,
    LOCATION_PLOT_LABEL_NEG_LOW,
    LOCATION_PLOT_LABEL_NEG_MID,
    LOCATION_PLOT_LABEL_NEG_HIGH,
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


def plot_count_clustermap():
    """
    Plots a clustermap of all samples using normalised counts.
    """
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
    relevant_peak_ids = []
    relevant_gene_symbols = []
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.startswith(
                "annotated_differential_accessibility"
            ) and file.endswith(".csv"):
                da_table_path = os.path.join(root, file)
                print(
                    f"\tProcessing {da_table_path}...",
                    flush=True,
                )
                da_table = pd.read_csv(
                    os.path.join(root, file),
                    sep=",",
                    header=0,
                    index_col=0,
                    encoding="utf-8",
                )
                filtered_da_table = da_table[da_table["padj"] <= 0.05]
                print(
                    f"\tFound {len(filtered_da_table.index)} significantly different regions. Selecting top peaks...",
                    flush=True,
                )
                table_of_relevant_peaks = filtered_da_table.sort_values(
                    "padj", ascending=True
                ).head(MAX_RELEVANT_GENES_PER_DA)
                relevant_peak_ids.extend(table_of_relevant_peaks.index)
                relevant_gene_symbols.extend(table_of_relevant_peaks["Gene Name"])
    relevant_peaks = OrderedDict(zip(relevant_peak_ids, relevant_gene_symbols))
    print(f"Selected {len(relevant_peaks)} relevant peaks...", flush=True)

    print("Plotting clustermap...", flush=True)
    # The set is sorted into a list afterwards to allow consistent assignment
    # colour assignment per label.
    unique_col_names = sorted(set(count_matrix_conditions), key=str.casefold)
    col_colour_palette = sns.color_palette(
        palette=COLOUR_PALETTE, n_colors=len(unique_col_names)
    )
    colour_mapping = OrderedDict(zip(unique_col_names, col_colour_palette))
    # Assigns the column names colour based on their respective condition.
    column_colours = {}
    for index_col, name_col in enumerate(count_matrix_col_names):
        column_colours[name_col] = colour_mapping[count_matrix_conditions[index_col]]

    plot = sns.clustermap(
        count_matrix.loc[list(relevant_peaks.keys())],
        cmap="vlag",
        col_colors=pd.Series(column_colours, name="Sample type"),
        cbar_kws={"label": "standardised rlog counts"},
        standard_scale=0,
        center=0.5,
    )
    plot.ax_row_dendrogram.set_visible(False)

    # Plots a legend for the per sample type / condition colour code.
    handles = [Patch(facecolor=colour_mapping[name]) for name in colour_mapping]
    plt.legend(
        handles,
        colour_mapping,
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

    print("Saving clustermap peak annotation...", flush=True)
    key_list = list(relevant_peaks.keys())
    value_list = list(relevant_peaks.values())
    reordered_dataframe = pd.DataFrame(
        data={
            "Peak ID": [key_list[i] for i in plot.dendrogram_row.reordered_ind],
            "Gene Name": [value_list[i] for i in plot.dendrogram_row.reordered_ind],
        }
    )
    reordered_dataframe.to_csv(
        os.path.join(
            MOUNT_PATHS["output"],
            "clustermap_peak_annotation.csv",
        ),
        sep=",",
        encoding="utf-8",
    )
    plt.close()


def annotation_to_plot_name(da_row):
    """
    Converts annotations to plotable categories.
    """
    annotation_name = da_row["Annotation"]
    for key in ANNOTATION_RENAMING_MAP:
        if (
            annotation_name is not None
            # Filters empty columns that are interpreted as NAs.
            and not isinstance(annotation_name, float)
            and annotation_name.startswith(key)
        ):
            return ANNOTATION_RENAMING_MAP[key]
    # Returns NULL if no matching annotation can be found.
    return None


def lfc_to_category(da_row):
    """
    Converts log2-fold-changes to plotable categories.
    """
    lfc = float(da_row["log2FoldChange"])
    for label, cutoff in LOCATION_PLOT_BINS.items():
        if lfc >= cutoff:
            return label
    # Returns NULL if no matching cutoff.
    return None


def sign_for_counts(cutoff_value) -> float:
    """
    Returns -1 for negative cutoff values and +1 for zero and positive cutoffs.
    """
    if np.sign(cutoff_value) < 0.0:
        return -1.0
    else:
        return 1.0


def plot_genomic_region_barplots():
    print("Processing differential accessibility data...", flush=True)
    for root, dirs, files in os.walk(INPUT_FOLDER):
        for file in files:
            if file.startswith(
                "annotated_differential_accessibility"
            ) and file.endswith(".csv"):
                da_file = os.path.join(root, file)
                print(f"Processing {da_file}...", flush=True)
                output_path = os.path.join(
                    MOUNT_PATHS["output"],
                    file.replace(
                        "annotated_differential_accessibility",
                        "genomic_regions_barplot",
                        1,
                    ).removesuffix(".csv")
                    + ".svg",
                )
                da_table = pd.read_csv(
                    os.path.join(da_file),
                    sep=",",
                    header=0,
                    index_col=0,
                    encoding="utf-8",
                )
                print("\tFiltering data...", flush=True)
                filtered_da_table = da_table[da_table["padj"] <= 0.05][
                    ["log2FoldChange", "Annotation"]
                ].copy()
                filtered_da_table["plot_annotation"] = filtered_da_table.apply(
                    annotation_to_plot_name, axis=1
                )
                filtered_da_table["lfc_category"] = filtered_da_table.apply(
                    lfc_to_category, axis=1
                )
                filtered_da_table = filtered_da_table[
                    (filtered_da_table["plot_annotation"].notnull())
                    & (filtered_da_table["lfc_category"].notnull())
                ]
                if len(filtered_da_table.index) == 0:
                    print(
                        "\tNo differentially accessible peaks found. Skipping sample...",
                        flush=True,
                    )
                    continue
                # Creates columns that contain the counts of
                # all the log2-fold-change cutoff categories.
                count_dictionary = {}
                for label, cutoff in LOCATION_PLOT_BINS.items():
                    count_dictionary[label] = filtered_da_table[
                        filtered_da_table["lfc_category"] == label
                    ].groupby(["plot_annotation"]).size() * sign_for_counts(cutoff)
                location_summaries = pd.DataFrame(count_dictionary)
                # Missing values are just zero counts.
                location_summaries.fillna(value=0, inplace=True)
                location_summaries = location_summaries.reindex(
                    columns=LOCATION_PLOT_LABEL_ORDER_STACK,
                    copy=False,
                )
                print("\tCreating genomic location plot...", flush=True)
                stacked_barplot_ax = location_summaries.plot.bar(
                    color=LOCATION_PLOT_COLOURS, rot=0, stacked=True, figsize=(9, 5)
                )
                # Hides the axis labels.
                stacked_barplot_ax.xaxis.label.set_visible(False)
                # Reorders legend labels.
                stacked_barplot_handles, stacked_barplot_labels = (
                    stacked_barplot_ax.get_legend_handles_labels()
                )
                stacked_barplot_handle_map = dict(
                    zip(stacked_barplot_labels, stacked_barplot_handles)
                )
                ordered_handles = list(
                    map(
                        stacked_barplot_handle_map.get, LOCATION_PLOT_LABEL_ORDER_LEGEND
                    )
                )
                stacked_barplot_ax.legend(
                    handles=ordered_handles, labels=LOCATION_PLOT_LABEL_ORDER_LEGEND
                )
                # Converts negative Y axis labels to positve counts.
                stacked_barplot_ax.set_yticklabels(
                    [f"{abs(x):0.0f}" for x in stacked_barplot_ax.get_yticks()]
                )
                # Saves the plot.
                stacked_barplot_figure = stacked_barplot_ax.get_figure()
                stacked_barplot_figure.savefig(output_path)
                plt.close(stacked_barplot_figure)


# Runs the visualisation functions.
plot_count_clustermap()
plot_genomic_region_barplots()
