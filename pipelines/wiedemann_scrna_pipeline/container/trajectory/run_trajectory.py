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
from py_monocle import (
    learn_graph,
    order_cells,
    compute_cell_states,
    regression_analysis,
    differential_expression_genes,
)
from scipy import sparse

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
KEY_ROOT_CELL = "iroot"
KEY_DIFFUSION_MAP_X = "X_diffmap"
KEY_PSEUDOTIME_DIFFUSION = "dpt_pseudotime"
KEY_PSEUDOTIME_MONOCLE = "monocle_pseudotime"
KEY_CELLSTATE_MONOCLE = "monocle_cellstate"

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


def key_cluster_pseudotime(obs_name, pseudotime_key):
    """
    Returns a cluster specific version of the pseudotime, key.
    """
    return f"{pseudotime_key}_{cluster_name_to_cluster(obs_name)}"


print("Loading cluster tree data...", flush=True)
adata = anndata.read_h5ad(os.path.join(INPUT_FOLDER, "cluster_relation.h5ad"))
print("Calculating diffusion maps...", flush=True)
sc.tl.diffmap(adata, n_comps=50)
adata.uns[KEY_ROOT_CELL] = adata.obsm[KEY_DIFFUSION_MAP_X][:, 3].argmin()
sc.tl.dpt(adata)
print("Plotting data...", flush=True)
fig = sc.pl.umap(
    adata,
    color=KEY_PSEUDOTIME_DIFFUSION,
    legend_loc="on data",
    show=False,
    return_fig=True,
)
fig.savefig(os.path.join(MOUNT_PATHS["output"], "diffusion_pseudotime_umap.svg"))
plt.close(fig)

print("Setting up Monocle...", flush=True)
monocle_umap_data = adata.obsm["X_umap"]
monocle_expression_matrix = sparse.csr_matrix(adata.X)

print("Analysing clusters...", flush=True)
for obs_name in adata.obs:
    if is_cluster_obs(obs_name=obs_name):
        cluster = cluster_name_to_cluster(obs_name=obs_name)
        print(f"\tProcessing {cluster} clusters...", flush=True)
        if int(cluster) < 3:
            print(f"\t\tNot enough clusters. Skipping...", flush=True)
        else:
            print("\t\tRunning Monocle...", flush=True)
            key_monocle_pseudotime = key_cluster_pseudotime(
                obs_name=obs_name, pseudotime_key=KEY_PSEUDOTIME_MONOCLE
            )
            key_monocle_cellstate = key_cluster_pseudotime(
                obs_name=obs_name, pseudotime_key=KEY_CELLSTATE_MONOCLE
            )
            projected_points, mst, centroids = learn_graph(
                matrix=monocle_umap_data, clusters=adata.obs[obs_name]
            )
            adata.obs[key_monocle_pseudotime] = order_cells(
                monocle_umap_data,
                centroids,
                mst=mst,
                projected_points=projected_points,
                root_cells=adata.uns[KEY_ROOT_CELL],
            )
            cell_states, _ = compute_cell_states(monocle_umap_data, centroids, mst)
            adata.obs[key_monocle_cellstate] = cell_states

            print("\t\tPlotting Monocle...", flush=True)

            fig_monocel_principal, ax_monocole_principal = plt.subplots(figsize=(8, 6))
            sc.pl.umap(
                adata,
                color=obs_name,
                legend_loc="on data",
                show=False,
                return_fig=False,
                ax=ax_monocole_principal,
            )
            edges = np.array(mst.nonzero()).T
            for edge in edges:
                sns.lineplot(
                    x=centroids[edge, 0],
                    y=centroids[edge, 1],
                    sort=False,
                    color="black",
                    linewidth=1,
                    ax=ax_monocole_principal,
                )
            ax_monocole_principal.set_title("Principal graph")
            fig_monocel_principal.savefig(
                os.path.join(MOUNT_PATHS["output"], f"monocle_principal_{cluster}.svg")
            )
            plt.close(fig_monocel_principal)

            fig_monocel_cellstate, ax_monocole_cellstate = plt.subplots(figsize=(8, 6))
            sc.pl.umap(
                adata,
                color=key_monocle_cellstate,
                legend_loc="on data",
                show=False,
                return_fig=False,
                ax=ax_monocole_cellstate,
            )
            edges = np.array(mst.nonzero()).T
            for edge in edges:
                sns.lineplot(
                    x=centroids[edge, 0],
                    y=centroids[edge, 1],
                    sort=False,
                    color="black",
                    linewidth=1,
                    ax=ax_monocole_cellstate,
                )
            ax_monocole_cellstate.set_title("Principal graph with ceLtb4r1ll states")
            fig_monocel_cellstate.savefig(
                os.path.join(MOUNT_PATHS["output"], f"monocle_cellstate_{cluster}.svg")
            )
            plt.close(fig_monocel_cellstate)

            fig_monocel_pseudotime = sc.pl.umap(
                adata,
                color=key_monocle_pseudotime,
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            fig_monocel_pseudotime.savefig(
                os.path.join(MOUNT_PATHS["output"], f"monocle_pseudotime_{cluster}.svg")
            )
            plt.close(fig_monocel_pseudotime)

            print("\t\tPlotting Monocle differential gene expression...", flush=True)
            morans_i_scores, pvalues, adjusted_pvalues = differential_expression_genes(
                monocle_expression_matrix,
                projected_cells=projected_points,
            )
            morani_results = pd.DataFrame(
                {
                    "Moran's I score": morans_i_scores,
                    "p-value": pvalues,
                    "Adjusted p-value": adjusted_pvalues,
                },
                index=pd.Series(adata.var_names, name="Gene"),
            ).sort_values("Moran's I score", ascending=False)
            morani_results.to_csv(
                os.path.join(
                    MOUNT_PATHS["output"], f"monocle_differential_gene_expression_{cluster}.csv"
                ),
                sep=",",
                encoding="utf-8",
            )
            fig_monocel_dge = sc.pl.umap(
                adata,
                color=morani_results.index[:32],
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            fig_monocel_dge.savefig(
                os.path.join(MOUNT_PATHS["output"], f"monocle_dge_{cluster}.svg")
            )
            plt.close(fig_monocel_dge)

            print("\t\tRunning PAGA...", flush=True)
            sc.tl.paga(adata, groups=obs_name)
            print("\t\tPlotting PAGA...", flush=True)
            ax_paga = sc.pl.paga(adata, show=False)
            fig_paga_graph = ax_paga.get_figure()
            fig_paga_graph.savefig(
                os.path.join(MOUNT_PATHS["output"], f"paga_graph_{cluster}.svg")
            )
            plt.close(fig_paga_graph)

            print("\t\tPlotting PAGA comparisons...", flush=True)
            with plt.rc_context({"figure.figsize": (8, 8)}):
                ax_paga = sc.pl.paga_compare(
                    adata, show=False, threshold=0.10, edge_width_scale=0.25
                )
                fig_paga_graph = ax_paga[0].get_figure()
                fig_paga_graph.savefig(
                    os.path.join(MOUNT_PATHS["output"], f"paga_compare_{cluster}.svg")
                )
                plt.close(fig_paga_graph)
