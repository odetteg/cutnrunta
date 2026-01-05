import os
import sys
from pathlib import Path
from constants.common import *


configfile: "config/config.yaml"


include: "workflows/rules/qc.smk"
include: "workflows/rules/idx.smk"
include: "workflows/rules/map.smk"
include: "workflows/rules/multiqc.smk"
include: "workflows/rules/plots.smk"
include: "workflows/rules/read_len.smk"
include: "workflows/rules/bigwig.smk"
include: "workflows/rules/deeptools.smk"
include: "workflows/rules/bam2bed.smk"
include: "workflows/rules/peaks.smk"
include: "workflows/rules/peaks_qc.smk"
include: "workflows/rules/bam_merge.smk"
include: "workflows/rules/bamcompare.smk"
include: "workflows/rules/h3k27_shift.smk"
include: "constants/common.smk"


rule all:
    input:
        targets(),
