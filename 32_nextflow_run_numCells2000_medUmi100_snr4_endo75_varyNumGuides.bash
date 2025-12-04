#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="numCells2000_medUmi100_snr4_endo75_varyNumGuides"

    source include/nextflow_run_common.sh

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
