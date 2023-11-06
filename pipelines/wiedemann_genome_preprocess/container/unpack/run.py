#!/usr/bin/python
"""This module unpacks the genome file."""

import json
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

# Extract file.
for root, dirs, files in os.walk(MOUNT_PATHS['input']):
    for file in files:
        if file.casefold().endswith(".fa.gz"):
            file_input_path = os.path.join(root, file)
            full_command = (f"gzip -dkc {file_input_path} "
            f"> {MOUNT_PATHS['output']}/genome.fa")
            exit_code = os.waitstatus_to_exitcode(os.system(full_command))
            if exit_code != 0:
                sys.exit(exit_code)
