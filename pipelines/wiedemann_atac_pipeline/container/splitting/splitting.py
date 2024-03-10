#!/usr/bin/python
"""This module splits fragments into nucleosome specific chunks."""

import json
import math
import multiprocessing
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["alignment_filtering"]

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

# Iterates over all sample directories and processes them conserving the directory structure.
BAM_SUFFIX = ".bam"
for root, dirs, files in os.walk(INPUT_FOLDER):
    if len(files) > 0:
        for file in files:
            if file.endswith(BAM_SUFFIX):
                file_base_name = file.removesuffix(BAM_SUFFIX)
                file_base_path = os.path.join(root, file_base_name)
                file_base_output_path = os.path.join(
                    MOUNT_PATHS["output"],
                    file_base_path.removeprefix(INPUT_FOLDER + "/"),
                )
                #   1 -  99 bp: nucleosome free region
                # 180 - 247 bp: mononucleosomal region
                # 315 - 473 bp: dinucleosomal region
                # 558 - 615 bp: trinucleosomal region
                full_commands = [
                    (
                        f"samtools view -@ {threads} -h {file_base_path}{BAM_SUFFIX} "
                        "| awk '{if (($9 < 100 && $9 > -100 && $9 != 0) || $1 ~ /^@/) {print $0}}' "
                        f"| samtools sort -@ {threads} "
                        f"-O bam -o {file_base_output_path}_nucelosomefree.bam -"
                    ),
                    (
                        f"samtools view -@ {threads} -h {file_base_path}{BAM_SUFFIX} "
                        "| awk '{if ((($9 <= 247 && $9 >= 180) || ($9 >= -247 && $9 <= -180)) || $1 ~ /^@/) {print $0}}' "
                        f"| samtools sort -@ {threads} "
                        f"-O bam -o {file_base_output_path}_mononucelosomal.bam -"
                    ),
                    (
                        f"samtools view -@ {threads} -h {file_base_path}{BAM_SUFFIX} "
                        "| awk '{if ((($9 <= 473 && $9 >= 315) || ($9 >= -473 && $9 <= -315)) || $1 ~ /^@/) {print $0}}' "
                        f"| samtools sort -@ {threads} "
                        f"-O bam -o {file_base_output_path}_dinucelosomal.bam -"
                    ),
                    (
                        f"samtools view -@ {threads} -h {file_base_path}{BAM_SUFFIX} "
                        "| awk '{if ((($9 <= 615 && $9 >= 558) || ($9 >= -615 && $9 <= -558)) || $1 ~ /^@/) {print $0}}' "
                        f"| samtools sort -@ {threads} "
                        f"-O bam -o {file_base_output_path}_trinucelosomal.bam -"
                    ),
                ]
                for full_command in full_commands:
                    print(f"Running command: {full_command}")
                    os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
                    exit_code = os.waitstatus_to_exitcode(os.system(full_command))
                    if exit_code != 0:
                        sys.exit(exit_code)

for root, dirs, files in os.walk(MOUNT_PATHS["output"]):
    if len(files) > 0:
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
