#!/usr/bin/python
"""This module removes blacklisted regions."""

import json
import math
import multiprocessing
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))
BLACKLIST_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "blacklist.bed")
BAM_SUFFIX = ".bam"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            print(f"Processing file {os.path.join(root, file)}...", flush=True)
            file_base_name = file.removesuffix(BAM_SUFFIX)
            file_base_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_path.removeprefix(INPUT_FOLDER + "/"),
            )
            os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)

            print("Removing blacklisted regions...", flush=True)
            subprocess.run(
                (
                    f"bedtools intersect -abam {file_base_path}{BAM_SUFFIX} -b {BLACKLIST_PATH} -wa -v "
                    f"| samtools sort -@ {threads} -O bam -o {file_base_output_path}{BAM_SUFFIX} -"
                ),
                shell=True,
                check=True,
            )
            subprocess.run(
                f"samtools index -@ {threads} {file_base_output_path}{BAM_SUFFIX}",
                shell=True,
                check=True,
            )
            subprocess.run(
                (
                    f"samtools flagstat -@ {threads} -O json {file_base_output_path}{BAM_SUFFIX} > "
                    f"{file_base_output_path}.flagstat"
                ),
                shell=True,
                check=True,
            )
