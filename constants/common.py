import os
import pandas as pd
from pathlib import Path
import yaml
from typing import List, Tuple


#Accessing the config file and setting the directories
BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"
with open(CONFIG_FILE_PATH, 'r') as yamfile:
    config = yaml.safe_load(yamfile)
RAW_DATA_DIR = config.get('RAW_DATA_DIR')
REF_DIR = config.get('REF_DIR')
SO_DIR = config.get('SO_DIR')
RESULTS_DIR = config.get('RESULTS_DIR')
readS = ["R1", "R2"]
fastqc_dir = str(RESULTS_DIR) + "/fastqc"
logs_dir = str(SO_DIR) + "/logs"


samples_df = pd.read_csv(os.path.join(BASE_DIR, "samples.txt"), header=None)
all_ids = []
all_meta_ids = []
for id_ in samples_df.iloc[:, 0].to_list():
    filename = Path(id_).name
    parts = filename.split("_", 1)
    sample_id = parts[0]
    meta_id = parts[1]
    all_ids.append(sample_id)
    all_meta_ids.append(meta_id)
sample_ids = list(set(all_ids))
meta_ids = list(set(all_meta_ids))
root_meta_ids = []

for meta_id in meta_ids:
    p = Path(meta_id)
    while p.suffix:
        p = p.with_suffix('')  
    root_meta_ids.append(p.name)

print(fastqc_dir)

