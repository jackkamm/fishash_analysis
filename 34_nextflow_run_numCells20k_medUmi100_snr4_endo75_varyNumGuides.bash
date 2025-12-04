#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="numCells20k_medUmi100_snr4_endo75_varyNumGuides"

    source include/nextflow_run_common.sh

    # note demuxEM requires a lot of memory in this scenario
    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --skipCrispatBig \
             --skipCleanser \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
