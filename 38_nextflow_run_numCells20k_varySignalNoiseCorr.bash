#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="numCells20k_varySignalNoiseCorr"

    source include/nextflow_run_common.sh

    # only using this run to compare Simpson's correction vs without,
    # so skip slower methods

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --maxMemFactor 2 \
             --skipCrispatPoisson \
             --skipCrispatNegBinom \
             --skipSceptre \
             --skipCleanser \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
