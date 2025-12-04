
if [ ! -d "outs" ]; then
    echo 'Please create directory (or symlink) named "outs/" to store outputs'
    exit 1
fi

source 00_setup_env.sh

sh 10_simulate.sh

sh 20_download_barnyard_data.sh
Rscript 21_process_barnyard_data.R

bash 30_assign_guides.bash

Rscript 41_plot_simulation_results.R
Rscript 42_plot_barnyard_results.R
