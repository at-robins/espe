#!/usr/bin/python
"""This module runs the trimming process."""

import json
import math
import multiprocessing
import os
import pysam
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

# Sets environment variables.
threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

remove_invalid_reads = os.environ.get("REMOVE_INVALID")
if remove_invalid_reads == "true":
    print("Activated filter: remove invalid reads", flush=True)

remove_mitochondrial_reads = os.environ.get("REMOVE_M")
if remove_mitochondrial_reads == "true":
    print("Activated filter: remove mitochondrial reads", flush=True)

quality_filtering = os.environ.get("QUALITY_FILTER")
if quality_filtering is not None:
    print(
        f"Activated filter: minimum mapping quality of {quality_filtering}", flush=True
    )


def passes_quality_filtering(bam_row):
    """
    Checks if quality filtering should be performed and if a specific alignment
    passes the filtering.
    """
    return quality_filtering is None or bam_row.mapping_quality >= int(
        quality_filtering
    )


def passes_invalid_reads(bam_row):
    """
    Checks if invalid alignments should be removed and if a specific alignment
    passes the filtering.
    """
    # Samtools flag -F 2828 -f 3
    return remove_invalid_reads != "true" or (
        (
            not bam_row.mate_is_unmapped
            and not bam_row.is_unmapped
            and not bam_row.is_supplementary
            and not bam_row.is_secondary
            and not bam_row.is_qcfail
        )
        and bam_row.is_paired
        and bam_row.is_proper_pair
    )


def passes_mito_reads(bam_row):
    """
    Checks if mitochondrial reads should be removed and if a specific alignment
    passes the filtering.
    """
    return (
        remove_mitochondrial_reads != "true" or bam_row.reference_name.lower() != "chrm"
    )


# Iterates over all sample directories and processes them conserving the directory structure.
INPUT_SUFFIX = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(INPUT_SUFFIX):
            file_base_name = file.removesuffix(INPUT_SUFFIX)
            file_base_input_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_input_path.removeprefix(INPUT_FOLDER + "/"),
            )
            os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
            print(f"Filtering file {file_base_input_path}{INPUT_SUFFIX}...", flush=True)
            with pysam.AlignmentFile(
                f"{file_base_input_path}{INPUT_SUFFIX}", mode="rb", threads=threads
            ) as bam_file_in:
                with pysam.AlignmentFile(
                    f"{file_base_output_path}{INPUT_SUFFIX}",
                    mode="wb",
                    header=bam_file_in.header,
                    threads=threads,
                ) as bam_file_out:
                    fetched_alignments = bam_file_in.fetch(until_eof=True)
                    for alignment_row in fetched_alignments:
                        if (
                            passes_invalid_reads(alignment_row)
                            and passes_mito_reads(alignment_row)
                            and passes_quality_filtering(alignment_row)
                        ):
                            bam_file_out.write(alignment_row)
