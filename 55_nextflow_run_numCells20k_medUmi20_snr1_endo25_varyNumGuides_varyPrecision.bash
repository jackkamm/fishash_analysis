#!/bin/bash
set -euxo pipefail
{
    SIM_NAME="numCells20k_medUmi20_snr1_endo25_varyNumGuides"
    RUN_NAME="$SIM_NAME"_varyPrecision

    source include/nextflow_run_common.sh

    mkdir -p $OUTS/results/$RUN_NAME

    nextflow ../../vary_precision.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$SIM_NAME \
             -profile $NXF_PROFILE \
             --cleanserMemFactor 8 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/simulations/$SIM_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
