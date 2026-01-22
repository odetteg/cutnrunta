import pandas as pd
import json
import subprocess
import os

bam = "/shared/projects/cutnrunta/so_/results/bam/mark_remove_dups/TBL3-H2K119Ub-2_S14_L001.deduped.bam"

cmd = f"samtools view -F 4 {bam} | head -n 1000000 | cut -f 10 | awk '{print length}' | sort | uniq -c | sort -nr | head -n 1"

import subprocess

res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
print(res.stdout.strip())  # should print: 8 59

mode_len = int(res.stdout.strip().split()[-1])
print(mode_len)  # prints: 59
