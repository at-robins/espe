#!/usr/bin/python
"""This module runs the differential accessibility analysis."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
COUNT_MATRIX_PATH = os.path.join(INPUT_FOLDER, "counts/count_matrix.txt")

print("\tLoading R dependencies...")
importr("DESeq2")
print("\tRunning R...")
edger_function = ro.r(
        """
        function(data, sample_reference, sample_test, output_path){


fc_pe <- featureCounts(file_list, annot.ext = regions_to_count, isPairedEnd = FALSE, nthreads = 12)
counts <- fc_pe$counts
colnames(counts) <- c(
  "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
  "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
  "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"
)

group <- factor(c(
  "control", "IL12/18", "IL2/15", "IL12/18/2/15", "TGFb", "IL12/18/TGFb", "IL2/15/TGFb", "IL12/18/2/15/TGFb",
  "control", "IL12/18", "IL2/15", "IL12/18/2/15", "TGFb", "IL12/18/TGFb", "IL2/15/TGFb", "IL12/18/2/15/TGFb",
  "control", "IL2/15", "IL12/18", "IL12/18/2/15", "TGFb", "IL2/15/TGFb", "IL12/18/2/15/TGFb", "IL12/18/TGFb"
))
meta_data <- data.frame(group, row.names = colnames(counts))
atac_dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta_data, design = ~group)
atac_dds <- DESeq(atac_dds, parallel = TRUE)
atac_rlog <- rlog(atac_dds)
# atac_rlog
svg(filename = "13_pca.svg")
plotPCA(atac_rlog, intgroup = "group", ntop = nrow(atac_rlog))
dev.off()

get_result <- function(dds, sample_condition, reference_condition) {
  res <- results(dds, contrast = c("group", sample_condition, reference_condition))
  res["condition"] <- sample_condition
  return(res)
}

result_ctl <- get_result(atac_dds, "TGFb", "control")
result_12_18 <- get_result(atac_dds, "IL12/18/TGFb", "IL12/18")
result_2_15 <- get_result(atac_dds, "IL2/15/TGFb", "IL2/15")
result_12_18_2_15 <- get_result(atac_dds, "IL12/18/2/15/TGFb", "IL12/18/2/15")

results_total <- rbind(
  result_ctl,
  result_12_18,
  result_2_15,
  result_12_18_2_15
)

row.names(results_total) <- c(
    paste(row.names(result_ctl), "TGFb", sep = custom_seperator),
    paste(row.names(result_12_18), "IL12/18/TGFb", sep = custom_seperator),
    paste(row.names(result_2_15), "IL2/15/TGFb", sep = custom_seperator),
    paste(row.names(result_12_18_2_15), "IL12/18/2/15/TGFb", sep = custom_seperator)
)

gene_name <- c()
chromosome <- c()
gene_id <- c()
peak_start <- c()
peak_end <- c()
i_gi <- 0
for (gene_info in strsplit(row.names(results_total), custom_seperator)) {
  # Throw an error if the seperation yields to many results.
  if (length(gene_info) != 6) {
    stop(gene_info)
  }
  # Tracks the current index
  i_gi <- i_gi + 1
  gene_name[i_gi] <- gene_info[1]
  chromosome[i_gi] <- gene_info[2]
  gene_id[i_gi] <- gene_info[3]
  peak_start[i_gi] <- gene_info[4]
  peak_end[i_gi] <- gene_info[5]
}
results_total$geneName <- gene_name
results_total$chromosome <- chromosome
results_total$geneId <- gene_id
results_total$peakStart <- peak_start
results_total$peakEnd <- peak_end


write.csv(
  results_total,
  file = "14_hits_unfiltered.csv"
)


           # Creates an edgeR object with counts and grouping factor.
            y <- DGEList(assay(data, "X"), group = colnames(data))
            # Filters out genes with low counts.
            cat("\\tFeatures and samples before subsetting: ", dim(y)[1], " x ", dim(y)[2], "\\n", sep="")
            keep <- filterByExpr(y)
            y <- y[keep, , keep.lib.sizes=FALSE]
            cat("\\tFeatures and samples after subsetting: ", dim(y)[1], " x ", dim(y)[2], "\\n", sep="")
            # Performs normalisation.
            y <- calcNormFactors(y)
            # Creates a vector that is concatentation of condition and variables that we will later use with contrasts
            sample <- colData(data)$sample
            # Creates the design matrix.
            design <- model.matrix(~ 0 + sample)
            # Estimates dispersion.
            y <- estimateDisp(y, design = design)
            # Fits the model.
            fit <- glmQLFit(y, design)
            cat("\\tPlotting data...\\n", sep="")
            svg(paste(output_path, "mds_plot.svg", sep = "/"))
            plotMDS(y)
            dev.off()
            svg(paste(output_path, "bcv_plot.svg", sep = "/"))
            plotBCV(y)
            dev.off()
            # eval workaround as make makeContrasts does not accept a string variable.
            contrast_string = paste("sample", sample_test, " - sample", sample_reference, sep = "")
            cmd <- paste("sample_contrast <- makeContrasts(", contrast_string, ", levels = y$design)", sep ='"')
            eval(parse(text = cmd))
            qlf <- glmQLFTest(fit, contrast=sample_contrast)
            # Returns all of the DE genes and calculates Benjamini-Hochberg adjusted FDR.
            tt <- topTags(qlf, n = Inf)
            tt <- tt$table
            write.csv(tt, paste(output_path, "differential_gene_expression.csv", sep = "/"), row.names=TRUE)
            svg(paste(output_path, "smear_plot.svg", sep = "/"))
            plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<=0.05)])
            dev.off()
        }
        """
    )
    edger_function(
        pseudobulk_adata,
        VALUE_PSEUDOBULK_SAMPLE_REFERENCE,
        VALUE_PSEUDOBULK_SAMPLE_TEST,
        output_path,
    )