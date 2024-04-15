#!/usr/bin/python
"""This module runs the trimming process."""

import json
import math
import multiprocessing
import os
import pysam
import sys
import time

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["alignment"]

# If a specific environment variable is set, appends the respective option.

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

view_options = f"-@ {threads} -h "

remove_invalid_reads = os.environ.get("REMOVE_INVALID")
if remove_invalid_reads is not None and remove_invalid_reads == "true":
    view_options += "-F 2828 "

quality_filtering = os.environ.get("QUALITY_FILTER")
if quality_filtering is not None:
    view_options += f"-q {quality_filtering} "

remove_mitochondrial_reads = os.environ.get("REMOVE_M")


def passes_quality_filtering(bam_row):
    return quality_filtering is None or bam_row.mapping_quality >= int(
        quality_filtering
    )


def passes_invalid_reads(bam_row):
    return remove_invalid_reads != "true" or (
        not bam_row.mate_is_unmapped
        and not bam_row.is_unmapped
        and not bam_row.is_supplementary
        and not bam_row.is_secondary
        and not bam_row.is_qcfail
    )


def passes_mito_reads(bam_row):
    return (
        remove_mitochondrial_reads != "true" or bam_row.reference_name.lower() != "chrm"
    )


# Iterates over all sample directories and processes them conserving the directory structure.
INPUT_SUFFIX = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    if len(files) > 0:
        for file in files:
            if file.endswith(INPUT_SUFFIX):
                file_base_name = file.removesuffix(INPUT_SUFFIX)
                file_base_input_path = os.path.join(root, file_base_name)
                file_base_output_path = os.path.join(
                    MOUNT_PATHS["output"],
                    file_base_input_path.removeprefix(INPUT_FOLDER + "/"),
                )
                start = time.time()
                full_command = (
                    f"samtools view {view_options}{file_base_input_path}{INPUT_SUFFIX} "
                )
                remove_mitochondrial_reads = os.environ.get("REMOVE_M")
                if (
                    remove_mitochondrial_reads is not None
                    and remove_mitochondrial_reads == "true"
                ):
                    full_command += "| awk '{if($3 != \"chrM\"){print $0}}' "
                full_command += (
                    f"| samtools view -h -b - "
                    f"> {file_base_output_path}.bam"
                )
                print(f"Running command: {full_command}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)
                command_line = time.time()
                print(command_line - start, flush=True)
                with pysam.AlignmentFile(
                    f"{file_base_input_path}{INPUT_SUFFIX}", mode="rb", threads=threads
                ) as bam_file_in:
                    with pysam.AlignmentFile(
                        f"{file_base_output_path}_test{INPUT_SUFFIX}",
                        mode="wb",
                        header=bam_file_in.header,
                        threads=threads,
                    ) as bam_file_out:
                        fetched = bam_file_in.fetch(until_eof=True)
                        for x in fetched:
                            if (
                                passes_invalid_reads(x)
                                and passes_mito_reads(x)
                                and passes_quality_filtering(x)
                            ):
                                bam_file_out.write(x)
                        internal = time.time()
                        print(internal - command_line, flush=True)
