#!/usr/bin/python
"""This module calculates differentially expressed genes."""

import anndata
import anndata2ri
import csv
import json
import logging
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import warnings

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["pseudobulk_generation"] + "/"
CELL_TYPE_KEY = "celltype"
SAMPLE_KEY = "sample"
REPLICATE_KEY = "replicate"

# Setup of rpy2.
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()
ro.r(
    '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))'
)


# Setup of scanpy.
sc.settings.verbosity = 2


def differential_gene_expression(
    pseudobulk_adata, sample_type_a, sample_type_b, cell_type, output_folder_path
):
    """
    Performes differential gene expression analysis.
    """
    print(
        f"Comparing {sample_type_a} and {sample_type_b} for cell type {cell_type}",
        flush=True,
    )
    os.makedirs(output_folder_path, exist_ok=True)

    # Subsetting data.
    sample_a_mask = adata.obs[SAMPLE_KEY] == sample_type_a
    sample_b_mask = adata.obs[SAMPLE_KEY] == sample_type_b
    cell_type_mask = adata.obs[CELL_TYPE_KEY] == cell_type
    data_mask = np.logical_and(
        np.logical_or(sample_a_mask, sample_b_mask), cell_type_mask
    )
    adata_subset = adata[data_mask]

    print("\tLoading R dependencies...")
    importr("Seurat")
    importr("edgeR")
    print("\tRunning R...")
    edger_function = ro.r(
        """
        function(data, sample_a, sample_b, output_path){
            # create an edgeR object with counts and grouping factor
            y <- DGEList(assay(data, "X"), group = colnames(data))
            # filter out genes with low counts
            print("Dimensions before subsetting:")
            print(dim(y))
            keep <- filterByExpr(y)
            y <- y[keep, , keep.lib.sizes=FALSE]
            print("Dimensions after subsetting:")
            print(dim(y))
            # normalize
            y <- calcNormFactors(y)
            # create a vector that is concatentation of condition and cell type that we will later use with contrasts
            replicate <- colData(data)$replicate
            sample <- colData(data)$sample
            # create a design matrix: here we have multiple donors so also consider that in the design matrix
            design <- model.matrix(~ 0 + sample)
            # estimate dispersion
            y <- estimateDisp(y, design = design)
            # fit the model
            fit <- glmQLFit(y, design)
            print("\tPlotting data...")
            svg(paste(output_path, "mds_plot.svg", sep = "/"))
            plotMDS(y)
            dev.off()
            svg(paste(output_path, "bcv_plot.svg", sep = "/"))
            plotBCV(y)
            dev.off()
            print(colnames(y$design))
            contrast_string = paste("sample", sample_a, "-sample", sample_b, sep = "")
            print(contrast_string)
            myContrast <- makeContrasts(contrast_string, levels = y$design)
            qlf <- glmQLFTest(fit, contrast=myContrast)
            # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
            tt <- topTags(qlf, n = Inf)
            tt <- tt$table
            write.csv(tt, paste(output_path, "differentiall_gene_expression.csv", sep = "/"), row.names=TRUE)
            return(list("fit"=fit, "design"=design, "y"=y))
        }
        """
    )
    edger_output = edger_function(
        adata_subset, sample_type_a, sample_type_b, output_folder_path
    )


print("Reading pseudobulk data...")
adata = anndata.read_h5ad(f"{INPUT_FOLDER}pseudobulk.h5ad")

# Replaces invalid characters in sample names.
adata.obs[SAMPLE_KEY] = list(
    map(lambda val: val.replace(" ", "_"), adata.obs[SAMPLE_KEY])
)
adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype("category")

# Runs differential gene expression analysis.
differential_gene_expression(
    adata,
    sample_type_a="control_tissue",
    sample_type_b="tumourous_tissue",
    cell_type="NK cell",
    output_folder_path=MOUNT_PATHS["output"],
)
