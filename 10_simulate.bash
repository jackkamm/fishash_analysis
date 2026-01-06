#!/bin/bash
set -euxo pipefail
{
    source include/env_vars.sh
    OUTS=$(realpath $OUTS)

    mkdir -p $OUTS/simulations

    ## small test scenario to make sure things run
    #Rscript simulate/simulate_test.R $OUTS/simulations/test

    ## 20k cells, high-gRNA regime, varying number of guides
    #Rscript \
    #    simulate/simulate_numCells20k_medUmi100_snr4_endo75_varyNumGuides.R \
    #    $OUTS/simulations/numCells20k_medUmi100_snr4_endo75_varyNumGuides

    ## 20k cells, low-gRNA regime, varying number of guides
    #Rscript \
    #    simulate/simulate_numCells20k_medUmi20_snr1_endo25_varyNumGuides.R \
    #    $OUTS/simulations/numCells20k_medUmi20_snr1_endo25_varyNumGuides

    # varying MOI
    Rscript \
        simulate/simulate_numCells20k_numGuides200_varyMOI.R \
        $OUTS/simulations/numCells20k_numGuides200_varyMOI

    ## varying signal-noise correlation
    #Rscript \
    #    simulate/simulate_numCells20k_varySignalNoiseCorr.R \
    #    $OUTS/simulations/numCells20k_varySignalNoiseCorr
}
