
First, create a subfolder (or symlink) named `outs/` in the current directory,
where the outputs will be stored.

Then, run each of the numbered scripts in turn, or run them all together using
`run_all.sh`.

Prerequisites:
* Nextflow
* Conda, configured to install from conda-forge and bioconda repos
* Slurm, needed to run the nextflow jobs
  * Alternatively, edit [env_vars.sh](env_vars.sh) and set
    `NX_PROFILE="local"` to configure Nextflow to run all jobs locally
