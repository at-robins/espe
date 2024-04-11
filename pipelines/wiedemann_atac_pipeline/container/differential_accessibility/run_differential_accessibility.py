#!/usr/bin/python
"""This module runs the differential accessibility analysis."""

import csv
import json
import os
import pathvalidate
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

from rpy2.robjects.packages import importr

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
COUNT_MATRIX_PATH = os.path.join(INPUT_FOLDER, "counts/count_matrix_final.txt")
SAMPLE_COMPARISON_PATH = os.path.join(MOUNT_PATHS["input"], "sample_comparison.csv")

print("Loading count matrix...", flush=True)
peak_ids = []
count_matrix_values_by_row = []
groups = []
with open(
    COUNT_MATRIX_PATH, newline="", mode="rt", encoding="utf-8"
) as count_matrix_in:
    count_matrix_reader = csv.reader(
        count_matrix_in, dialect="unix", delimiter="\t", quotechar='"'
    )
    for i_matrix, row_matrix in enumerate(count_matrix_reader):
        peak_row = row_matrix[0]
        sample_rows = row_matrix[1 : len(row_matrix)]
        if i_matrix == 0:
            # Reads the condition names.
            groups = sample_rows
        else:
            # Reads the condition values.
            peak_ids.append(peak_row)
            for sample_value in sample_rows:
                count_matrix_values_by_row.append(int(sample_value))

print("Parsing information for sample comparisons...")
samples_test = []
samples_reference = []
sample_comparisons_output_suffixes = []
with open(SAMPLE_COMPARISON_PATH, newline="", encoding="utf-8") as csvfile:
    info_reader = csv.DictReader(csvfile, dialect="unix", delimiter=",", quotechar='"')
    for row in info_reader:
        sample_reference = row["reference sample"]
        sample_test = row["test sample"]
        samples_reference.append(sample_reference)
        samples_test.append(sample_test)
        sample_comparisons_output_suffixes.append(
            pathvalidate.sanitize_filename(f"{sample_test}__vs__{sample_reference}")
        )

r_count_matrix = ro.r.matrix(
    ro.IntVector(count_matrix_values_by_row),
    nrow=len(peak_ids),
    ncol=len(groups),
    byrow=True,
)

print("Loading R dependencies...", flush=True)
importr("DESeq2")
print("Running R...", flush=True)
deseq_function = ro.r(
    """
    function(
        data_matrix,
        peak_ids,
        groups,
        conditions_test,
        conditions_reference,
        output_path,
        output_suffixes
    ) {
        cat("Preprocessing data matrix...\\n", sep="")
        # Names rows and columns of the data matrix.
        rownames(data_matrix) <- peak_ids
        colnames(data_matrix) <- 1:ncol(data_matrix)
        meta_data <- data.frame(factor(groups), row.names = colnames(data_matrix))
        colnames(meta_data) <- c("groups")
        atac_dds <- DESeqDataSetFromMatrix(
            countData = data_matrix,
            colData = meta_data,
            design = ~groups
        )
        atac_dds <- DESeq(atac_dds, parallel = TRUE)

        cat("Regularising count matrix...\\n", sep="")
        atac_rlog <- rlog(atac_dds)

        cat("Plotting PCA...\\n", sep="")
        svg(filename = "PCA.svg")
        plotPCA(atac_rlog, intgroup = "groups", ntop = nrow(atac_rlog))
        dev.off()

        cat("Performing differential accessiblity analysis...\\n", sep="")
        for (i in 1:length(conditions_test)) {
            cat(
                "\\tComparing sample ",
                conditions_test[i],
                " to reference ",
                conditions_reference[i],
                "...\\n",
                sep=""
            )
            atac_result <- results(
                atac_dds,
                contrast = c("groups", conditions_test[i], conditions_reference[i])
            )
            atac_result["conditionSample"] <- conditions_test[i]
            atac_result["conditionReference"] <- conditions_reference[i]
            write.csv(
                atac_result,
                file = paste(
                    output_path,
                    "/",
                    "differential_accessibility_",
                    output_suffixes[i],
                    ".csv",
                    sep = ""
                ),
                row.names=TRUE
            )
        }
        cat("Done.\\n", sep="")
    }
    """
)
deseq_function(
    r_count_matrix,
    peak_ids,
    groups,
    samples_test,
    samples_reference,
    MOUNT_PATHS["output"],
    sample_comparisons_output_suffixes,
)
