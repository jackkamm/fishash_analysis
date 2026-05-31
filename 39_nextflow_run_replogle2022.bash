#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="replogle2022"

    source include/nextflow_run_common.sh

    # Note: Nearly all the crispat-negbinom jobs finished in 24h, but
    # a small handful (out of 273) did not. Shame to lose
    # crispat-negbinom for such a relatively small number of failures,
    # so I added option --slurmClusterOptionsVeryLong and reran it
    # with 3day queue instead of the usual 24 hours

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --skipCombineAssignments \
             --slurmClusterOptionsVeryLong '--nodes=1 --qos=3d' \
             --crispatBigMemFactor 2 \
             --cleanserMemFactor 4 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/replogle2022/processed/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
