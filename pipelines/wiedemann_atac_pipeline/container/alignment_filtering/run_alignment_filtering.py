#!/usr/bin/python
"""This module runs the trimming process."""

import json
import math
import multiprocessing
import os
import pysam

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
    input_files = [file for file in files if file.endswith(INPUT_SUFFIX)]
    if len(input_files) > 0:
        output_directory = os.path.join(
            MOUNT_PATHS["output"],
            os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
        )
        os.makedirs(output_directory, exist_ok=True)
        output_file_paths = []

        for file in input_files:
            input_file_path = os.path.join(root, file)
            output_file_path = os.path.join(output_directory, file)
            output_file_paths.append(output_file_path)
            print(f"Filtering file {input_file_path}...", flush=True)
            with pysam.AlignmentFile(
                f"{input_file_path}", mode="rb", threads=threads
            ) as bam_file_in:
                with pysam.AlignmentFile(
                    f"{output_file_path}",
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

        # Normally there should not be multiple sequencing files
        # per sample, but in rare cases samples might be sequenced
        # on different lanes etc. and need to be merged for
        # subsequent analysis.
        if len(output_file_paths) > 1:
            print(f"\tMerging {output_file_paths}...", flush=True)
            merged_file_path = os.path.join(output_directory, f"merged{INPUT_SUFFIX}")
            pysam.merge("@", threads, "-o", merged_file_path, *output_file_paths)
            print("\tDeleting merged files...", flush=True)
            for output_file in output_file_paths:
                os.remove(output_file)
            print("\tSorting merged file...", flush=True)
            merged_sorted_file_path = os.path.join(
                output_directory, f"merged_sorted{INPUT_SUFFIX}"
            )
            pysam.sort(
                "-O",
                "bam",
                "@",
                threads,
                "-o",
                merged_sorted_file_path,
                merged_file_path,
            )
            print("\tDeleting unsorted file...", flush=True)
            os.remove(merged_file_path)
