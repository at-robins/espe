#!/usr/bin/python
"""This module performs pathway analysis."""

import anndata
import csv
print("Loading decoupler...", flush=True)
import decoupler
import json
import numpy as np
import os
import pandas as pd
import pathvalidate
import scanpy as sc
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
KEY_PATHWAY_DB_GENESET = "source"
KEY_PATHWAY_DB_GENESYMBOL = "target"
FILTER_PATHWAY_DB_MIN = 5
FILTER_PATHWAY_DB_MAX = 1500
NETWORK_PATH_HUMAN_PROGENY = "/pathway_networks/human.pkl"
NETWORK_PATH_MOUSE_PROGENY = "/pathway_networks/mouse.pkl"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.set_figure_params(
    scanpy=True,
    # In case of bitmap exports, use high quality.
    dpi=300,
    dpi_save=300,
    # Export as SVG.
    format="svg",
    vector_friendly=True,
    # Use transparent background.
    transparent=True,
    facecolor=None,
    # Remove frames.
    frameon=False,
)
sc.settings.figdir = MOUNT_PATHS["output"]


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


def perform_enrichment_analysis(
    file_path_input, pathway_db, output_folder, suffix
):
    """
    Runs the ULM.
    """
    os.makedirs(output_folder, exist_ok=True)
    print(f"Processing file {file_path_input}", flush=True)
    print("\tReading data...", flush=True)
    adata = anndata.read_h5ad(file_path_input)
    print("\tRunning ULM...", flush=True)
    decoupler.mt.ulm(data=adata, net=pathway_db, tmin=3)
    score = decoupler.pp.get_obsm(adata, key="score_ulm")

    print("\tPlotting data...", flush=True)
    all_umap_samples = score.var_names
    chunk_size = 4 * 8
    for chunk_index in range(0, len(all_umap_samples), chunk_size):
        fig = sc.pl.umap(
            score,
            color=all_umap_samples[chunk_index : chunk_index + chunk_size],
            wspace=1,
            show=False,
            return_fig=True,
            cmap="RdBu_r",
            vcenter=0,
        )
        fig.savefig(os.path.join(output_folder_path, f"pathways_{suffix}_{chunk_index}.svg"))
        plt.close(fig)



print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    pathway_database_reactome = parse_gmt(GMT_PATH_HUMAN_REACTOME)
    pathway_database_go = parse_gmt(GMT_PATH_HUMAN_GO)
    pathway_database_hallmark = parse_gmt(GMT_PATH_HUMAN_HALLMARK)
    pathway_database_progeny = pd.read_pickle(NETWORK_PATH_HUMAN_PROGENY)
else:
    print("\tUsing mouse data...", flush=True)
    pathway_database_reactome = parse_gmt(GMT_PATH_MOUSE_REACTOME)
    pathway_database_go = parse_gmt(GMT_PATH_MOUSE_GO)
    pathway_database_hallmark = parse_gmt(GMT_PATH_MOUSE_HALLMARK)
    pathway_database_progeny = pd.read_pickle(NETWORK_PATH_MOUSE_PROGENY)

# Iterates over all directories and loads sample information.
print("Searching for data...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.casefold().endswith(".h5ad"):
            file_path_input = os.path.join(root, file)
            output_folder_path = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            os.makedirs(output_folder_path, exist_ok=True)

            print("Using PROGENY database...", flush=True)
            perform_enrichment_analysis(
                file_path_input=file_path_input,
                pathway_db=pathway_database_progeny,
                output_folder=output_folder_path,
                suffix="progeny",
            )
            print("Using hallmark database...", flush=True)
            perform_enrichment_analysis(
                file_path_input=file_path_input,
                pathway_db=pathway_database_hallmark,
                output_folder=output_folder_path,
                suffix="hallmark",
            )
            print("Using reactome database...", flush=True)
            perform_enrichment_analysis(
                file_path_input=file_path_input,
                pathway_db=pathway_database_reactome,
                output_folder=output_folder_path,
                suffix="reactome",
            )
            print("Using gene ontology database...", flush=True)
            perform_enrichment_analysis(
                file_path_input=file_path_input,
                pathway_db=pathway_database_go,
                output_folder=output_folder_path,
                suffix="go",
            )
