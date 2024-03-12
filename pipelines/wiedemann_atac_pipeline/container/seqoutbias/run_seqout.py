#!/usr/bin/python
"""This module corrects transposase induced bias."""

import json
import math
import multiprocessing
import os
import subprocess
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = MOUNT_PATHS["dependencies"]["splitting"]
STRAND_MASK_PLUS = "NXNXXXCXXNNXNNNXXN"
STRAND_MASK_MINUS = "NXXNNNXNNXXCXXXNXN"
GENOME_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "genome.fa.gz")
READ_SIZE = os.environ.get("READ_LENGTH")
TALLYMER_PATH = os.path.join(
    MOUNT_PATHS["dependencies"]["tallymer"], f"genome.tal_{READ_SIZE}.gtTxt.gz"
)
SEQTABLE_PATH_PLUS = os.path.join(
    MOUNT_PATHS["dependencies"]["seqtable"], "seqtable_plus.tbl"
)
SEQTABLE_PATH_MINUS = os.path.join(
    MOUNT_PATHS["dependencies"]["seqtable"], "seqtable_minus.tbl"
)
BAM_SUFFIX = ".bam"
BW_SUFFIX = ".bigWig"

threads = math.floor(multiprocessing.cpu_count() * 0.8)
if threads < 1:
    threads = 1

# Iterates over all sample directories and processes them conserving the directory structure.
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.endswith(BAM_SUFFIX):
            file_base_name = file.removesuffix(BAM_SUFFIX)
            file_base_path = os.path.join(root, file_base_name)
            file_base_output_path = os.path.join(
                MOUNT_PATHS["output"],
                file_base_path.removeprefix(INPUT_FOLDER + "/"),
            )
            bam_path_plus = f"{file_base_output_path}_plus{BAM_SUFFIX}"
            bam_path_minus = f"{file_base_output_path}_minus{BAM_SUFFIX}"
            bw_path_plus = f"{file_base_output_path}_plus{BW_SUFFIX}"
            bw_path_minus = f"{file_base_output_path}_minus{BW_SUFFIX}"
            full_commands = [
                f"samtools view -@ {threads} -bh -F 16 {file_base_path}{BAM_SUFFIX} -o {bam_path_plus}",
                f"samtools view -@ {threads} -bh -f 16 {file_base_path}{BAM_SUFFIX} -o {bam_path_minus}",
                (
                    f"seqOutBias {GENOME_PATH} {bam_path_plus} "
                    f"--bw={bw_path_plus} --out={SEQTABLE_PATH_PLUS} "
                    f"--skip-bed --shift-counts --kmer-mask {STRAND_MASK_PLUS} "
                    f"--read-size={READ_SIZE} --tallymer={TALLYMER_PATH}"
                ),
                (
                    f"seqOutBias {GENOME_PATH} {bam_path_minus} "
                    f"--bw={bw_path_minus} --out={SEQTABLE_PATH_MINUS} "
                    f"--skip-bed --shift-counts --kmer-mask {STRAND_MASK_MINUS} "
                    f"--read-size={READ_SIZE} --tallymer={TALLYMER_PATH}"
                ),
                f"bigWigMerge {bw_path_plus} {bw_path_minus} {file_base_output_path}.bedGraph",
                # Cleans up temporary data.
                f"rm {bam_path_plus} {bam_path_minus}",
            ]
            for full_command in full_commands:
                print(f"Running command: {full_command}", flush=True)
                os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
                subprocess.run(
                    full_command,
                    shell=True,
                    check=True,
                )
