#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="replogle2022"

    source include/nextflow_run_common.sh

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --skipCrispatNegBinom \
             --crispatBigMemFactor 2 \
             --cleanserMemFactor 4 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/replogle2022/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
