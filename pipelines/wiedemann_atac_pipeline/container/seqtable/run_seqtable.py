#!/usr/bin/python
"""This module creates files required for seqOutBias."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
STRAND_MASK_PLUS = "NXNXXXCXXNNXNNNXXN"
STRAND_MASK_MINUS = "NXXNNNXNNXXCXXXNXN"
GENOME_PATH = os.path.join(MOUNT_PATHS["globals"]["GENOME"], "genome.fa.gz")
READ_SIZE = os.environ.get("READ_LENGTH")
TALLYMER_PATH = os.path.join(
    MOUNT_PATHS["dependencies"]["tallymer"], f"genome.tal_{READ_SIZE}.gtTxt.gz"
)

subprocess.run(
    (
        "seqOutBias seqtable "
        f"{os.path.join(MOUNT_PATHS['globals']['GENOME'], 'genome.fa.gz')} "
        f"--out={os.path.join(MOUNT_PATHS['output'],'seqtable_plus.tbl')} "
        f"--kmer-mask={STRAND_MASK_PLUS} "
        f"--read-size={READ_SIZE} --tallymer={TALLYMER_PATH}"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)

subprocess.run(
    (
        "seqOutBias seqtable "
        f"{os.path.join(MOUNT_PATHS['globals']['GENOME'], 'genome.fa.gz')} "
        f"--out={os.path.join(MOUNT_PATHS['output'],'seqtable_minus.tbl')} "
        f"--kmer-mask={STRAND_MASK_MINUS} "
        f"--read-size={READ_SIZE} --tallymer={TALLYMER_PATH}"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
