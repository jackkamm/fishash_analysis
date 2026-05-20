#!/usr/bin/env python

import os
import re
import sys

with open("include/env_vars.sh") as f:
    exec(f.read())

DATA_DIR = os.path.join(OUTS, "replogle2022")
RAW_DIR = os.path.join(DATA_DIR, "raw/K562_gwps_other")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")

if not os.path.exists(RAW_DIR):
    sys.exit(f"Please download files to {RAW_DIR} from https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_MTX_files/20127869/1")

mtx_re = re.compile(r'(.*)_matrix.mtx.gz')

batch_names = []
for fname in os.listdir(RAW_DIR):
    matched = mtx_re.match(fname)
    if matched:
        batch_names.append(matched.group(1))

if not batch_names:
    sys.exit(f"Did not find any mtx files in {RAW_DIR}, please make sure the data files are extracted there")

# HACK: 273 is hardcoded into 23_process_replogle2022.slurm as the
# number of batches, so I assert that it's the length here
assert len(batch_names) == 273

os.makedirs(PROCESSED_DIR, exist_ok=True)

with open(os.path.join(PROCESSED_DIR, "sample_sheet.csv"), "w") as f:
    print("sim_label,path", file=f)
    for b in batch_names:
        rds = os.path.join(PROCESSED_DIR, b + ".Rds")
        rds = os.path.realpath(rds)
        print(f"{b},{rds}", file=f)
