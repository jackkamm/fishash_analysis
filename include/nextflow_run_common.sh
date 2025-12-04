# preamble to be sourced in all nextflow_run scripts
# before sourcing this, need to set the RUN_NAME env var

source include/env_vars.sh
OUTS=$(realpath $OUTS)
NXF_WORKDIR=$(realpath $NXF_WORKDIR)
mkdir -p $OUTS/results

# run in a separate subfolder so we can run it in parallel with
# other nextflow jobs
mkdir -p nextflow_run/$RUN_NAME
cd nextflow_run/$RUN_NAME
