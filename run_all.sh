
# Recommended: 2 cpus, 64GB memory
./10_simulate.bash

# Recommended: 2 cpus, 4GB memory
./20_download_barnyard_data.bash
# Recommended: 8 cpus, 32GB memory
Rscript 21_process_barnyard_data.R

# Recommended: 2 cpus, 16GB memory.
#
# Running these in sequence may take a long time (days).  Consider
# running each of these in a separate screen/tmux pane in parallel
./31_nextflow_run_test.bash
./34_nextflow_run_numCells20k_medUmi100_snr4_endo75_varyNumGuides.bash
./35_nextflow_run_numCells20k_medUmi20_snr1_endo25_varyNumGuides.bash
./36_nextflow_run_numCells20k_numGuides200_varyMOI.bash
./37_nextflow_run_cleanser_barnyard_data.bash
./38_nextflow_run_numCells20k_varySignalNoiseCorr.bash

# Recommended: 4 cpus, 16GB memory.
Rscript 61_plot_simulation_results.R
Rscript 62_plot_barnyard_results.R
