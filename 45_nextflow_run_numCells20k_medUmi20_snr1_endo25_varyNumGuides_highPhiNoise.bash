#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="numCells20k_medUmi20_snr1_endo25_varyNumGuides_highPhiNoise"

    source include/nextflow_run_common.sh

    mkdir -p $OUTS/results/$RUN_NAME

    # 20-20k guide scenarios

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --crispatBigMemFactor 2 \
             --cleanserMemFactor 4 \
             --demuxemMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet_leq20k.csv \
             --outdir $OUTS/results/$RUN_NAME/leq20k

    # big 80k guide scenario
    #
    # skip crispat big models due to their long running time.

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --slurmClusterOptionsVeryLong '--nodes=1 --qos=3d' \
             --cleanserMemFactor 4 \
             --demuxemMemFactor 4 \
             --crispatBigMemFactor 4 \
             --maxMemFactor 2 \
             --sample_sheet $OUTS/simulations/$RUN_NAME/sample_sheet_80k.csv \
             --outdir $OUTS/results/$RUN_NAME/80k
}
