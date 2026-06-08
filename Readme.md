
This repo is to reproduce the results in our [fishash manuscript](https://www.biorxiv.org/content/10.64898/2026.01.22.701179).

The [master](https://github.com/jackkamm/fishash_analysis/tree/master) branch of this repo is for the original (v1) preprint results. The updated results for the revision are in the [revision](https://github.com/jackkamm/fishash_analysis/tree/revision) branch.

Prerequisites:
* Nextflow
* Conda, configured to install from conda-forge and bioconda repos
* Slurm, needed to run the nextflow jobs
  * Alternatively, edit [include/env_vars.sh](include/env_vars.sh) and set
    `NXF_PROFILE="local"` to configure Nextflow to run all jobs locally

After satisfying the prerequisites, do the following:

0. Create a folder to store results, and make a symlink to it named
   `outs` in the current directory.
   ```
   mkdir /path/to/outs
   ln -s /path/to/outs outs
   ```
1. Create a folder to store intermediate outputs from nextflow,
   and make a symlink to it named `work` in the current directory.
   ```
   mkdir /path/to/work
   ln -s /path/to/work work
   ```
2. Create the conda environment `fishash_analysis` by running:
   ```
   source setup_env.sh
   ```
3. Activate the environment:
   ```
   conda activate fishash_analysis
   ```
4. If using slurm, you may need to edit the parameters
   `slurmClusterOptionsDefault` and `slurmClusterOptionsLong` in
   `nextflow.config` to set the options appropriate for your HPC Slurm
   system (such as queue, qos, etc).
5. Run each of the numbered scripts in turn.
   * Alternatively, run them all in a single script with `run_all.sh`
   * On a shared HPC environment, it may be preferable to manually run
     each script in turn (rather than directly using `run_all.sh`),
     due to different resource requirements and running times of each
     of the steps. Refer to `run_all.sh` on how to run each individual
     step, and see the comments for the recommended resources (cpus,
     memory) to allocate each step.
