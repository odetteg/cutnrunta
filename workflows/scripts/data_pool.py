import os
import sys
import polars as pl

print("Python executable:", sys.executable)


print("Polars version:", pl.__version__)


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

    df = pl.DataFrame(all_data)

    df.write_csv(agg_out)


input_file = sys.argv[1]
output_file = sys.argv[2]

aggregate_samtools_stats(input_file, output_file)
