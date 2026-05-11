#!/bin/bash
set -euxo pipefail
{
    SIM_NAME="numCells20k_medUmi100_snr4_endo75_varyNumGuides"
    RUN_NAME="$SIM_NAME"_varyPrecision

    source include/nextflow_run_common.sh

    mkdir -p $OUTS/results/$RUN_NAME

    # first do the 20-20k guide scenario

    nextflow ../../vary_precision.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$SIM_NAME \
             -profile $NXF_PROFILE \
             --cleanserMemFactor 4 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/simulations/$SIM_NAME/sample_sheet_leq20k.csv \
             --outdir $OUTS/results/$RUN_NAME/leq20k

    # if that successfully finishes, then do the big 80k guide scenario

    nextflow ../../vary_precision.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$SIM_NAME \
             -profile $NXF_PROFILE \
             --cleanserMemFactor 8 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/simulations/$SIM_NAME/sample_sheet_80k.csv \
             --outdir $OUTS/results/$RUN_NAME/80k
}
