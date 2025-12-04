#!/bin/bash
set -euxo pipefail
{
    source include/env_vars.sh
    FISHASH_ANALYSIS=$(pwd -P)

    mkdir -p $OUTS/barnyard_data/raw
    cd $OUTS/barnyard_data/raw
    cat $FISHASH_ANALYSIS/data/barnyard_data_urls.txt | xargs -I% wget %
}
