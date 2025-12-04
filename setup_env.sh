conda create -n fishash_analysis
conda activate fishash_analysis

# python>=3.11: needed for cleanser
# uv: needed for geomux
# splatter: for sceptre workaround
# cmdstanpy, muon: specific versions needed for cleanser
conda install python=3.11 r \
      uv \
      r-devtools \
      r-tidyverse \
      r-extradistr \
      bioconductor-summarizedexperiment \
      bioconductor-sparsematrixstats \
      r-argparse \
      r-reticulate \
      bioconductor-splatter \
      r-patchwork \
      r-pals \
      bioconductor-multiassayexperiment \
      bioconductor-singlecellexperiment \
      cmdstanpy=1.2 muon=0.1.7

# FIXME: Rather than installing uv within the conda install, it's
# probably better to just install uv outside the env and then install
# geomux thru that. Or to use pip to install geomux within the conda env?

pip install pegasusio demuxem

# modified version of crispat that comments out some of the plotting
# functionality which can occasionally cause crashes
pip install git+https://github.com/jackkamm/crispat.git@nohist

pip install git+https://github.com/Gersbachlab-Bioinformatics/CLEANSER.git

# run cleanser once to ensure models are precompiled and avoid race conditions
# https://github.com/Gersbachlab-Bioinformatics/CLEANSER/issues/10
cleanser -i data/cleanser_test.mtx -p 1 -c 1 --dc
cleanser -i data/cleanser_test.mtx -p 1 -c 1 --cs

Rscript -e 'devtools::install("~/devel/fishash", upgrade=FALSE)'
Rscript -e 'devtools::install_github("katsevich-lab/sceptre", upgrade=FALSE)'

uv tool install geomux

conda env export > environment.yaml

conda deactivate
