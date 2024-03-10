#!/usr/bin/python
"""This module creates files required for seqOutBias."""

import json
import os
import subprocess

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

subprocess.run(
    (
        "seqOutBias tallymer "
        f"{os.path.join(MOUNT_PATHS['dependencies']['get_genome'], 'genome.fa')} "
        "36"
    ),
    cwd=MOUNT_PATHS["output"],
    shell=True,
    check=True,
)
