#!/usr/bin/python
"""This module performs pseudotime and trajectory analysis."""

import anndata
import csv
import json
import numpy as np
import os
import pandas as pd
import pathvalidate
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
KEY_ROOT_CELL = "iroot"
KEY_DIFFUSION_MAP_X = "X_diffmap"
KEY_PSEUDOTIME_DIFFUSION = "dpt_pseudotime"
KEY_PSEUDOTIME_DIFFUSION_SCALED = "dpt_pseudotime_scaled"

# Setup of scanpy.
sc.settings.verbosity = 2
sc.settings.set_figure_params(
    scanpy=True,
    # # In case of bitmap exports, use high quality.
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
sc.settings.figdir = ""


def is_cluster_obs(obs_name):
    """
    Returns True if the observation name contains cluster information.
    """
    return obs_name.startswith("number_of_clusters_")


def cluster_name_to_cluster(obs_name):
    """
    Returns the cluster number from the observation name.
    """
    return obs_name.removeprefix("number_of_clusters_")


print("Loading cluster tree data...", flush=True)
adata = anndata.read_h5ad(os.path.join(INPUT_FOLDER, "cluster_relation.h5ad"))
print("Calculating diffusion maps...", flush=True)
sc.tl.diffmap(adata, n_comps=50)
adata.uns[KEY_ROOT_CELL] = adata.obsm[KEY_DIFFUSION_MAP_X][:, 3].argmin()
sc.tl.dpt(adata)
print("Plotting data...", flush=True)
# pdt_mean = np.mean(adata.obs[KEY_PSEUDOTIME_DIFFUSION])
# pdt_std = np.std(adata.obs[KEY_PSEUDOTIME_DIFFUSION])
# pdt_min = pdt_mean - 0.5 * pdt_std
# pdt_max = pdt_mean + 0.5 * pdt_std
# adata = adata[adata.obs[KEY_PSEUDOTIME_DIFFUSION] >= pdt_min]
# adata = adata[adata.obs[KEY_PSEUDOTIME_DIFFUSION] <= pdt_max]
# adata.obs[KEY_PSEUDOTIME_DIFFUSION_SCALED] = (adata.obs[KEY_PSEUDOTIME_DIFFUSION] - pdt_mean) / pdt_std
fig = sc.pl.umap(
    adata,
    color=KEY_PSEUDOTIME_DIFFUSION,
    legend_loc="on data",
    show=False,
    return_fig=True,
)
fig.savefig(os.path.join(MOUNT_PATHS["output"], "pseudotime_umap.svg"))
plt.close(fig)

for obs_name in adata.obs:
    if is_cluster_obs(obs_name=obs_name):
        cluster = cluster_name_to_cluster(obs_name=obs_name)
        print(f"\tProcessing {cluster} clusters...", flush=True)
        if int(cluster) < 3:
            print(f"\t\tNot enough clusters. Skipping...", flush=True)
        else:
            sc.tl.paga(adata, groups=obs_name)
            ax_paga = sc.pl.paga(adata, show=False)
            fig_paga_graph = ax_paga.get_figure()
            fig_paga_graph.savefig(
                os.path.join(MOUNT_PATHS["output"], f"paga_graph_{cluster}.svg")
            )
            plt.close(fig_paga_graph)

            ax_paga = sc.pl.paga_compare(adata, show=False)
            fig_paga_graph = ax_paga[0].get_figure()
            fig_paga_graph.savefig(
                os.path.join(MOUNT_PATHS["output"], f"paga_compare_{cluster}.svg")
            )
            plt.close(fig_paga_graph)
            # sc.tl.umap(
            #     adata,
            #     init_pos="paga",
            # )
            # fig_paga_umap = sc.pl.umap(
            #     adata,
            #     color=[KEY_PSEUDOTIME_DIFFUSION, obs_name, "Gzma"],
            #     legend_loc="on data",
            #     show=False,
            #     return_fig=True,
            # )
            # fig_paga_umap.savefig(
            #     os.path.join(MOUNT_PATHS["output"], f"paga_umap_{cluster}.svg")
            # )
            # plt.close(fig_paga_umap)

# print("\tWriting merged data to file...")
# output_file_path = os.path.join(MOUNT_PATHS["output"], "merged.h5ad")
# os.makedirs(MOUNT_PATHS["output"], exist_ok=True)
# adata_merged.write(output_file_path, compression="gzip")
