
sh 10_simulate.sh

sh 20_download_barnyard_data.sh
Rscript 21_process_barnyard_data.R

bash 30_assign_guides.bash

Rscript 41_plot_simulation_results.R
Rscript 42_plot_barnyard_results.R
