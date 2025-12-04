
Prerequisites:
* Nextflow
* Conda, configured to install from conda-forge and bioconda repos
* Slurm, needed to run the nextflow jobs
  * Alternatively, edit [env_vars.sh](env_vars.sh) and set
    `NX_PROFILE="local"` to configure Nextflow to run all jobs locally

After satisfying the prerequisites, do the following:

1. Create a subfolder (or symlink) named `outs/` in the current directory,
   where the outputs will be stored.
2. Create the conda environment `fishash_analysis` by running:
   ```
   source setup_env.sh
   ```
3. Activate the environment:
   ```
   conda activate fishash_analysis
   ```
4. Then run each of the numbered scripts in turn, or run them all using
   ```
   sh run_all.sh
   ```
