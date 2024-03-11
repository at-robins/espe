#!/usr/bin/python
"""This module splits fragments into nucleosome specific chunks."""

import json
import math
import multiprocessing
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["alignment_filtering"]
BASE_COMMAND = "java -jar /picard.jar"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

# Iterates over all sample directories and processes them conserving the directory structure.
BAM_SUFFIX = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            file_base_name = file.removesuffix(BAM_SUFFIX)
            file_base_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_path.removeprefix(INPUT_FOLDER + "/"),
            )
            full_commands = [
                (
                    f"{BASE_COMMAND} CollectInsertSizeMetrics --METRIC_ACCUMULATION_LEVEL ALL_READS --ASSUME_SORTED true "
                    f"--INPUT {file_base_path}{BAM_SUFFIX} "
                    f"--OUTPUT {file_base_output_path}_insert_size_metrics.txt "
                    f"--Histogram_FILE {file_base_output_path}_insert_size_metrics_histogram.pdf"
                ),
                (
                    f"{BASE_COMMAND} MarkDuplicates --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES true "
                    f"--INPUT {file_base_path}{BAM_SUFFIX} "
                    f"--OUTPUT {file_base_output_path}{BAM_SUFFIX} "
                    f"--METRICS_FILE {file_base_output_path}_metrics.txt"
                ),
            ]
            for full_command in full_commands:
                print(f"Running command: {full_command}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)

for root, dirs, files in os.walk(MOUNT_PATHS["output"]):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            file_base_name = file.removesuffix(BAM_SUFFIX)
            file_base_output_path = os.path.join(root, file_base_name)
            full_commands = [
                f"samtools index -@ {threads} {file_base_output_path}{BAM_SUFFIX}",
                (
                    f"samtools flagstat -@ {threads} -O json {file_base_output_path}{BAM_SUFFIX} > "
                    f"{file_base_output_path}.flagstat"
                ),
            ]
            for full_command in full_commands:
                print(f"Running command: {full_command}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)
