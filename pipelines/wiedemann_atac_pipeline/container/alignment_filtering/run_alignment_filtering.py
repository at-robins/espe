#!/usr/bin/python
"""This module runs the trimming process."""

import json
import math
import multiprocessing
import os
import sys

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
                    file_base_input_path.removeprefix(INPUT_FOLDER + "/")
                )
                full_command = f"samtools view {view_options}{file_base_input_path}{INPUT_SUFFIX} "
                remove_mitochondrial_reads = os.environ.get("REMOVE_M")
                if remove_mitochondrial_reads is not None and remove_mitochondrial_reads == "true":
                    full_command += "| awk '{if($3 != \"chrM\"){print $0}}' "
                full_command += (f"| samtools sort -@ {threads} "
                f"-O bam -o {file_base_output_path}.bam -")
                print(f"Running command: {full_command}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok = True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)
