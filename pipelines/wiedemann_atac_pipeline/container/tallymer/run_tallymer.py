#!/usr/bin/python
"""This module creates files required for seqOutBias."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))
READ_SIZE = os.environ.get("READ_LENGTH")

subprocess.run(
    (
        "seqOutBias tallymer "
        f"{os.path.join(MOUNT_PATHS['globals']['GENOME'], 'genome.fa.gz')} "
        f"{READ_SIZE}"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
