#!/bin/bash
set -euxo pipefail
{
    SIM_NAME="test"
    RUN_NAME="test_varyPrecision"

    source include/nextflow_run_common.sh

    nextflow ../../vary_precision.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --sample_sheet $OUTS/simulations/$SIM_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
