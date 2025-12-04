#!/bin/bash
set -euxo pipefail
{
    RUN_NAME="barnyard_data"

    source include/nextflow_run_common.sh

    # NOTE: A fake ground-truth was created so that the simulation
    # pipeline would run to completion (as the nextflow pipeline computes
    # the accuracy at the end assuming a known ground truth). In
    # particular, all entries with >0 counts were labeled as being truly
    # present. These accuracy metrics should be ignored as it is just a
    # hack to run the pipeline!  Instead, extract the full results and
    # manually compare them to the true species after running this (this
    # is done in 42_plot_barnyard_results.R).
    #
    # (TODO: refactor the pipeline so it can be run without requiring
    # ground truth and computing accuracy metrics)

    nextflow ../../main.nf \
             -resume \
             -work-dir $NXF_WORKDIR/$RUN_NAME \
             -profile $NXF_PROFILE \
             --sample_sheet $OUTS/barnyard_data/processed/sample_sheet.csv \
             --outdir $OUTS/results/$RUN_NAME
}
