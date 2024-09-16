#!/usr/bin/python
"""This module runs the peak annotation."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

# ANNOTATION_TTS = "TTS"
# ANNOTATION_PROMOTOR = "promotor"
# ANNOTATION_INTRON = "intron"
# ANNOTATION_EXON = "exon"
# ANNOTATION_INTERGENIC = "intergenic"
# ANNOTATION_NON_CODING = "non-coding"
# ANNOTATION_UTR_3 = "3\' UTR"
# ANNOTATION_UTR_5 = "5\' UTR"
# ANNOTATION_UNKNOWN = "unknown"

print("Loading pathway database...")
organism_env = os.environ.get("GLOBAL_ORGANISM")
if organism_env is not None and organism_env == "human":
    print("\tUsing human data...", flush=True)
    genome_id = "hg38"
else:
    print("\tUsing mouse data...", flush=True)
    genome_id = "mm10"

# Iterates over all sample directories and processes them conserving the directory structure.
print("Annotating peaks...", flush=True)
for root, dirs, files in os.walk(INPUT_FOLDER):
    for file in files:
        if file.startswith("merged."):
            file_path_input = os.path.join(root, file)
            directory_path_output = os.path.join(
                MOUNT_PATHS["output"],
                os.path.normpath(os.path.relpath(root, INPUT_FOLDER)),
            )
            file_path_output = os.path.join(directory_path_output, file)
            print(f"\tProcessing file {file_path_input}...", flush=True)
            os.makedirs(directory_path_output, exist_ok=True)

            subprocess.run(
                (
                    f"annotatePeaks.pl {file_path_input} "
                    f" {genome_id} > {file_path_output}"
                ),
                cwd=directory_path_output,
                shell=True,
                check=True,
            )
