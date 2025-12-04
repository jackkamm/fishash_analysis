source env_vars.sh

# small test scenario to make sure things run
Rscript simulate/simulate_test.R $GUIDEBENDER_SIM_DIR/test

# 2k cells, high-gRNA regime, varying number of guides
Rscript \
    simulate/simulate_numCells2000_medUmi100_snr4_endo75_varyNumGuides.R \
    $GUIDEBENDER_SIM_DIR/numCells2000_medUmi100_snr4_endo75_varyNumGuides

# 2k cells, low-gRNA regime, varying number of guides
Rscript \
    simulate/simulate_numCells2000_medUmi20_snr1_endo25_varyNumGuides.R \
    $GUIDEBENDER_SIM_DIR/numCells2000_medUmi20_snr1_endo25_varyNumGuides

# 20k cells, high-gRNA regime, varying number of guides
Rscript \
    simulate/simulate_numCells20k_medUmi100_snr4_endo75_varyNumGuides.R \
    $GUIDEBENDER_SIM_DIR/numCells20k_medUmi100_snr4_endo75_varyNumGuides

# 20k cells, low-gRNA regime, varying number of guides
Rscript \
    simulate/simulate_numCells20k_medUmi20_snr1_endo25_varyNumGuides.R \
    $GUIDEBENDER_SIM_DIR/numCells20k_medUmi20_snr1_endo25_varyNumGuides

# varying MOI
Rscript \
    simulate/simulate_numCells2000_numGuides100_varyMOI.R \
    $GUIDEBENDER_SIM_DIR/numCells2000_numGuides100_varyMOI
