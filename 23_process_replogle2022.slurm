#!/bin/bash
#SBATCH -J process_replogle2022
#SBATCH --qos=3h
#SBATCH --array=1-273
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -o outs/replogle2022/log/%x_%a.out
#SBATCH -e outs/replogle2022/log/%x_%a.err
#SBATCH --open-mode truncate
#SBATCH --mail-type=BEGIN,END,FAIL # Type of email notification- BEGIN,END,FAIL,ALL

source include/env_vars.sh

bin/make_replogle2022_summarized_experiment.R \
    --mtx_folder $OUTS/replogle2022/raw/K562_gwps_other \
    --sample_sheet_csv $OUTS/replogle2022/processed/sample_sheet.csv \
    --idx $SLURM_ARRAY_TASK_ID \
    --crispr_lib_csv data/k562_genomewide_library.csv.gz
