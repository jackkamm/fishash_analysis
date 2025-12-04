#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="numCells20k_medUmi20_snr1_endo25_varyNumGuides"

    source include/nextflow_run_common.sh

    # NOTE: do NOT use -resume, as we use this run for timing.
    nextflow ../../main.nf \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --maxMemFactor 2 \
             --skipCrispatBig \
             --skipCleanser \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
