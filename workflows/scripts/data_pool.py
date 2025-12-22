import os
from pathlib import Path
import pandas as pd
import sys


def aggregate_samtools_stats(samtools_txt, agg_out):
    with open(samtools_txt, "r") as f:
        samtools_files = [line.strip() for line in f if line.strip()]

    all_data = []
    for f_path in samtools_files:
        sample_name = os.path.basename(f_path).split(".")[0]
        with open(f_path, "r") as f:
            for line in f:
                if line.startswith("IS"):
                    parts = line.strip().split("\t")
                    all_data.append(
                        {
                            "sample": sample_name,
                            "size": int(parts[1]),
                            "count": int(parts[2]),
                            "in_pairs": int(parts[3]),
                            "out_pairs": int(parts[4]),
                            "other_pairs": int(parts[5]),
                        }
                    )

    df = pd.DataFrame(all_data)
    df.to_csv(agg_out, sep=",", index=False)


if __name__ == "__main__":
    samtools_txt = sys.argv[1]
    agg_out = sys.argv[2]

    aggregate_samtools_stats(samtools_txt, agg_out)
