
# Recommended: 2 cpus, 64GB memory
./10_simulate.bash

# Recommended: 2 cpus, 4GB memory
./20_download_barnyard_data.bash
# Recommended: 8 cpus, 32GB memory
Rscript 21_process_barnyard_data.R

# This part requires slurm to be setup
# TODO: Make it work for non-slurm systems
./22_prep_replogle2022_processing.py
sbatch < 23_process_replogle2022.slurm 

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

./39_nextflow_run_replogle2022.bash

./44_nextflow_run_numCells20k_medUmi100_snr4_endo75_varyNumGuides_highPhiNoise.bash
./45_nextflow_run_numCells20k_medUmi20_snr1_endo25_varyNumGuides_highPhiNoise.bash
./46_nextflow_run_numCells20k_numGuides200_varyMOI_highPhiNoise.bash

./54_nextflow_run_numCells20k_medUmi100_snr4_endo75_varyNumGuides_varyPrecision.bash
./55_nextflow_run_numCells20k_medUmi20_snr1_endo25_varyNumGuides_varyPrecision.bash
./56_nextflow_run_numCells20k_numGuides200_varyMOI_varyPrecision.bash

# Recommended: 4 cpus, 16GB memory.
Rscript 81_plot_simulation_results.R
Rscript 82_plot_barnyard_results.R
Rscript 83_plot_replogle2022_results.R
Rscript 84_plot_vary_precision.R
