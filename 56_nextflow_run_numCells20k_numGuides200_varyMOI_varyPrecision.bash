#!/bin/bash
set -euxo pipefail
{
    SIM_NAME="numCells20k_numGuides200_varyMOI"
    RUN_NAME="$SIM_NAME"_varyPrecision

    source include/nextflow_run_common.sh

    nextflow ../../vary_precision.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --maxMemFactor 2 \
             --crispatBigMemFactor 2 \
             --cleanserMemFactor 4 \
             --sample_sheet $OUTS/simulations/$SIM_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
