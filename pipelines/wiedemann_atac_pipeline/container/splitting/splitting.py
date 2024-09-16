#!/usr/bin/python
"""This module splits fragments into nucleosome specific chunks."""

import json
import math
import multiprocessing
import os
import pysam
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
INPUT_FOLDER = next(iter(MOUNT_PATHS["dependencies"].values()))

#   1 -  99 bp: nucleosome free region
REGION_NUCLEOSOME_FREE_MIN = 1
REGION_NUCLEOSOME_FREE_MAX = 99
# 180 - 247 bp: mononucleosomal region
REGION_MONONUCLEOSOMAL_MIN = 180
REGION_MONONUCLEOSOMAL_MAX = 247
# 315 - 473 bp: dinucleosomal region
REGION_DINUCLEOSOMAL_MIN = 315
REGION_DINUCLEOSOMAL_MAX = 473
# 558 - 615 bp: trinucleosomal region
REGION_TRINUCLEOSOMAL_MIN = 558
REGION_TRINUCLEOSOMAL_MAX = 615

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
            os.makedirs(os.path.dirname(file_base_output_path), exist_ok=True)
            print(f"Splitting file {file_base_path}{BAM_SUFFIX}...", flush=True)
            with pysam.AlignmentFile(
                f"{file_base_path}{BAM_SUFFIX}", mode="rb", threads=threads
            ) as bam_file_in:
                with pysam.AlignmentFile(
                    f"{file_base_output_path}_nucleosomefree{BAM_SUFFIX}",
                    mode="wb",
                    header=bam_file_in.header,
                    threads=threads,
                ) as bam_file_out_nucleosome_free, pysam.AlignmentFile(
                    f"{file_base_output_path}_mononucleosomal{BAM_SUFFIX}",
                    mode="wb",
                    header=bam_file_in.header,
                    threads=threads,
                ) as bam_file_out_mononucleosomal, pysam.AlignmentFile(
                    f"{file_base_output_path}_dinucleosomal{BAM_SUFFIX}",
                    mode="wb",
                    header=bam_file_in.header,
                    threads=threads,
                ) as bam_file_out_dinucleosomal, pysam.AlignmentFile(
                    f"{file_base_output_path}_trinucleosomal{BAM_SUFFIX}",
                    mode="wb",
                    header=bam_file_in.header,
                    threads=threads,
                ) as bam_file_out_trinucleosomal:
                    alignments_nucleosome_free = 0
                    alignments_mononucleosomal = 0
                    alignments_dinucleosomal = 0
                    alignments_trinucleosomal = 0
                    alignments_other = 0
                    fetched_alignments = bam_file_in.fetch(until_eof=True)
                    for alignment_row in fetched_alignments:
                        template_length = abs(alignment_row.template_length)
                        if (
                            template_length >= REGION_NUCLEOSOME_FREE_MIN
                            and template_length <= REGION_NUCLEOSOME_FREE_MAX
                        ):
                            bam_file_out_nucleosome_free.write(alignment_row)
                            alignments_nucleosome_free += 1
                        elif (
                            template_length >= REGION_MONONUCLEOSOMAL_MIN
                            and template_length <= REGION_MONONUCLEOSOMAL_MAX
                        ):
                            bam_file_out_mononucleosomal.write(alignment_row)
                            alignments_mononucleosomal += 1
                        elif (
                            template_length >= REGION_DINUCLEOSOMAL_MIN
                            and template_length <= REGION_DINUCLEOSOMAL_MAX
                        ):
                            bam_file_out_dinucleosomal.write(alignment_row)
                            alignments_dinucleosomal += 1
                        elif (
                            template_length >= REGION_TRINUCLEOSOMAL_MIN
                            and template_length <= REGION_TRINUCLEOSOMAL_MAX
                        ):
                            bam_file_out_trinucleosomal.write(alignment_row)
                            alignments_trinucleosomal += 1
                        else:
                            alignments_other += 1
                    print(f"\tNucleosome free alignments: {alignments_nucleosome_free}")
                    print(f"\tMononucleosomal alignments: {alignments_mononucleosomal}")
                    print(f"\tDinucleosomal alignments: {alignments_dinucleosomal}")
                    print(f"\tTrinucleosomal alignments: {alignments_trinucleosomal}")
                    print(f"\tOther alignments: {alignments_other}")

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
