#!/usr/bin/python
"""This module runs the samtools indexing."""

import json
import os
import sys

MOUNT_PATHS = json.loads(os.environ.get("MOUNT_PATHS"))

samtools_command = ("samtools faidx "
                    f"{MOUNT_PATHS['dependencies']['unpack']}/genome.fa "
                    f"-o {MOUNT_PATHS['output']}/genome.fa.fai")

exit_code = os.waitstatus_to_exitcode(os.system(samtools_command))
if exit_code != 0:
    sys.exit(exit_code)

cut_command = ("cut -f 1,2 "
                f"{MOUNT_PATHS['output']}/genome.fa.fai "
                f"> {MOUNT_PATHS['output']}/genome.chrom.sizes")

exit_code = os.waitstatus_to_exitcode(os.system(cut_command))
if exit_code != 0:
    sys.exit(exit_code)
