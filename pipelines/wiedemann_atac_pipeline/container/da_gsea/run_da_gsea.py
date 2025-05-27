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
GMT_PATH_HUMAN_REACTOME = "/gmts/reactome_human.gmt"
GMT_PATH_MOUSE_REACTOME = "/gmts/reactome_mouse.gmt"
GMT_PATH_HUMAN_GO = "/gmts/go_human.gmt"
GMT_PATH_MOUSE_GO = "/gmts/go_mouse.gmt"
GMT_PATH_HUMAN_HALLMARK = "/gmts/hallmark_human.gmt"
GMT_PATH_MOUSE_HALLMARK = "/gmts/hallmark_mouse.gmt"
KEY_PATHWAY_DB_GENESET = "geneset"
KEY_PATHWAY_DB_GENESYMBOL = "genesymbol"
FILTER_PATHWAY_DB_MIN = 15
FILTER_PATHWAY_DB_MAX = 500
PLOT_PATHWAYS_MAX = 30
FILTER_VALUE_ABSOLUTE = "abs"
KEY_SCORE = "scores"


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


def parse_da_csv(dge_csv_path, filter=None) -> pd.DataFrame:
    """
    Parses a CSV file containing differential gene expression data.
    """
    # Only keeps relevant columns.
    annotated_dars = pd.read_csv(
        dge_csv_path, sep=",", header=0, index_col=None, encoding="utf-8"
    )[["stat", "padj", "log2FoldChange", "Gene Name"]].copy()
    annotated_dars.sort_values("padj", ascending=True, inplace=True)
    total_peaks = len(annotated_dars.index)

    # Optionally ignore regulation.
    if filter == FILTER_VALUE_ABSOLUTE:
        annotated_dars["stat"] = np.abs(annotated_dars["stat"])

    # Only keeps relevant columns and renames the remaining ones.
    annotated_dars.drop(labels=["log2FoldChange", "padj"], axis=1, inplace=True)
    annotated_dars.rename(columns={"stat": KEY_SCORE}, inplace=True)

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


def perform_enrichment_analysis(
    dge_path, pathway_db, output_folder, suffix, filter=None
):
    """
    Performs the pathway enrichment analysis.
    """
    os.makedirs(output_folder, exist_ok=True)

    print("\tParsing data data...")
    cluster_genes = parse_da_csv(dge_csv_path=dge_path, filter=filter)

    print("\tRunning GSEA...", flush=True)
    scores, norm, pvals = decoupler.run_gsea(
        cluster_genes.T,
        pathway_db,
        source=KEY_PATHWAY_DB_GENESET,
        target=KEY_PATHWAY_DB_GENESYMBOL,
    )

    print("\tSaving results...", flush=True)
    gsea_results = (
        pd.concat({"score": scores.T, "norm": norm.T, "pval": pvals.T}, axis=1)
        .droplevel(level=1, axis=1)
        .assign(
            **{"abs_norm": lambda x: np.abs(x["norm"])},
            **{
                "regulation": lambda x: pd.Categorical(
                    map(
                        lambda z: ("less accessible" if z < 0.0 else "more accessible"),
                        np.sign(x["score"]),
                    ),
                    categories=["more accessible", "less accessible"],
                    ordered=True,
                )
            },
        )
        .sort_values(["pval", "abs_norm"])
    )

    if filter == FILTER_VALUE_ABSOLUTE:
        # Remove regulation column if using abolute ordering as it does not have any meaning.
        gsea_results.drop(labels=["regulation"], axis=1, inplace=False).to_csv(
            os.path.join(
                output_folder,
                pathvalidate.sanitize_filename(f"gsea_{suffix}_table.csv"),
            ),
            sep=",",
            encoding="utf-8",
        )
    else:
        gsea_results.to_csv(
            os.path.join(
                output_folder,
                pathvalidate.sanitize_filename(f"gsea_{suffix}_table.csv"),
            ),
            sep=",",
            encoding="utf-8",
        )

    print("\tWriting per gene data...", flush=True)
    gsea_results_filtered_sets = (
        gsea_results[gsea_results["pval"] <= 0.05]
        .sort_values("norm", ascending=False)
        .index.values
    )
    dge_lookup = parse_da_csv_for_extension(dge_csv_path=dge_path)
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

    print("\tPlotting data...", flush=True)
    gsea_results_filtered = gsea_results[gsea_results["pval"] <= 0.05]
    if filter == FILTER_VALUE_ABSOLUTE:
        # Sorted by signed NES, so only the top of the list is selected.
        gsea_results_filtered = (
            gsea_results_filtered[gsea_results_filtered["norm"] >= 0.0]
            .sort_values("norm", ascending=False)
            .head(PLOT_PATHWAYS_MAX)
        )
    else:
        # Sorted by absolute NES, so more and less accessible regions are selected separately.
        gsea_results_filtered = gsea_results_filtered.head(PLOT_PATHWAYS_MAX).sort_values(
            "norm", ascending=False
        )

    if len(gsea_results_filtered) == 0:
        print("\tNo significantly altered pathways. Skipping sample...", flush=True)
    else:
        height = (7 * (len(gsea_results_filtered) / PLOT_PATHWAYS_MAX)) + 1
        width = (max(map(len, gsea_results_filtered.index.values)) / 10) + 4
        fig, ax = plt.subplots(figsize=(width, height))
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
        if filter == FILTER_VALUE_ABSOLUTE:
            # There is no split by regulation when using absolute sorting,
            # so the legend should be removed.
            legend.remove()
        else:
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
    pathway_database_hallmark = parse_gmt(GMT_PATH_HUMAN_HALLMARK)
else:
    print("\tUsing mouse data...", flush=True)
    pathway_database_reactome = parse_gmt(GMT_PATH_MOUSE_REACTOME)
    pathway_database_go = parse_gmt(GMT_PATH_MOUSE_GO)
    pathway_database_hallmark = parse_gmt(GMT_PATH_MOUSE_HALLMARK)

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
            print("\tUsing hallmark database...", flush=True)
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_hallmark,
                output_folder=output_folder_path,
                suffix="hallmark_normal",
            )
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_hallmark,
                output_folder=output_folder_path,
                suffix="hallmark_ignore_regulation",
                filter=FILTER_VALUE_ABSOLUTE,
            )
            print("\tUsing reactome database...", flush=True)
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_reactome,
                output_folder=output_folder_path,
                suffix="reactome_normal",
            )
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_reactome,
                output_folder=output_folder_path,
                suffix="reactome_ignore_regulation",
                filter=FILTER_VALUE_ABSOLUTE,
            )
            print("\tUsing gene ontology database...", flush=True)
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_go,
                output_folder=output_folder_path,
                suffix="go_normal",
            )
            perform_enrichment_analysis(
                dge_path=dge_file,
                pathway_db=pathway_database_go,
                output_folder=output_folder_path,
                suffix="go_ignore_regulation",
                filter=FILTER_VALUE_ABSOLUTE,
            )
