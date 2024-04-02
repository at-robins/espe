#!/usr/bin/python
"""This module runs the alignment process."""

import json
import math
import multiprocessing
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["trimming"]

# The max fragment length of 615 bp allows to detect trinucleosomal regions.
options_bowtie = ("--end-to-end --no-mixed --dovetail --very-sensitive -X 615 "
f"-x {MOUNT_PATHS['globals']['GENOME']}/genome "
f"--met 10")

options_sort = ""
threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads > 0:
    options_bowtie += f" --threads {threads}"
    options_sort += f" -@ {threads}"


print("Specified bowtie2 options:" + options_bowtie)
print("Specified samtools options:" + options_sort)

# Iterates over all sample directories and processes them conserving the directory structure.
INPUT_SUFFIX_FORWARD = "_1_paired.fq.gz"
INPUT_SUFFIX_REVERSE = "_2_paired.fq.gz"
for root, dirs, files in os.walk(INPUT_FOLDER):
    if len(files) > 0:
        for file in files:
            if file.endswith(INPUT_SUFFIX_FORWARD):
                file_base_name = file.removesuffix(INPUT_SUFFIX_FORWARD)
                file_base_input_path = os.path.join(root, file_base_name)
                file_base_output_path = os.path.join(
                    MOUNT_PATHS["output"],
                    file_base_input_path.removeprefix(INPUT_FOLDER + "/")
                )
                full_command = (f"bowtie2 {options_bowtie} "
                f"--met-file {file_base_output_path}_metrics.txt "
                f"-1 {file_base_input_path}{INPUT_SUFFIX_FORWARD} "
                f"-2 {file_base_input_path}{INPUT_SUFFIX_REVERSE} | "
                f"samtools sort {options_sort} "
                f"-O bam -o {file_base_output_path}.bam -")
                print(f"Running command: {full_command}")
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok = True)
                exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                if exit_code != 0:
                    sys.exit(exit_code)
