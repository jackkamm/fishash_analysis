# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis pipeline for the fishash manuscript (https://www.biorxiv.org/content/10.64898/2026.01.22.701179v1). Fishash is a method for demultiplexing single-cell CRISPR experiments using guide RNA barcodes. This repo benchmarks fishash against competing methods on simulated and real data.

## Environment Setup

```bash
conda activate fishash_analysis
```

The full environment is specified in `environment.yaml`. To recreate from scratch, follow `setup_env.sh` (installs R 4.4, Python 3.11, and custom forks of crispat and CLEANSER).

The fishash R package itself lives at `~/devel/fishash` (local install, not on CRAN/Bioconductor).

## Running the Pipeline

The pipeline uses **Nextflow** with a SLURM backend. Each numbered bash script is a pipeline step:

```bash
# 1. Generate simulated data (2 CPUs, 64GB RAM)
./10_simulate.bash

# 2. Download and process real data
./20_download_barnyard_data.bash
Rscript 21_process_barnyard_data.R   # 8 CPUs, 32GB RAM

# 3. Run benchmarks (2 CPUs, 16GB RAM each — run in parallel in separate panes)
./31_nextflow_run_test.bash          # quick test
./34_nextflow_run_numCells20k_medUmi100_snr4_endo75_varyNumGuides.bash
./35_nextflow_run_numCells20k_medUmi20_snr1_endo25_varyNumGuides.bash
./36_nextflow_run_numCells20k_numGuides200_varyMOI.bash
./37_nextflow_run_cleanser_barnyard_data.bash
./38_nextflow_run_numCells20k_varySignalNoiseCorr.bash
./39_nextflow_run_replogle2022.bash

# 4. Generate plots (4 CPUs, 16GB RAM)
Rscript 61_plot_simulation_results.R
Rscript 62_plot_barnyard_results.R
Rscript 63_plot_replogle2022_results.R
```

Each `3x_nextflow_run_*.bash` script sources `include/nextflow_run_common.sh` and runs Nextflow in its own `nextflow_run/<RUN_NAME>/` subdirectory so multiple runs can execute in parallel.

Always pass `-resume` to Nextflow to restart from the last checkpoint (crispat's dask workers can timeout intermittently, requiring retries).

## Local vs SLURM Execution

Execution profile is controlled by `NXF_PROFILE` in `include/env_vars.sh` (default: `slurm`). To run locally, set `NXF_PROFILE="local"` there.

Output and work directories default to symlinks `outs/` and `work/`.

## Architecture

**Data flow:**
1. `bin/simulate_*.R` scripts generate `SummarizedExperiment` RDS files
2. `bin/convert_to_anndata.R` and `bin/convert_to_mtx.R` convert RDS → h5ad / MTX for Python tools
3. `main.nf` runs all demultiplexing methods in parallel
4. `bin/process_assignments.R` computes confusion matrices per method/replicate
5. `bin/combine_assignments.R` and `bin/combine_confusion.R` aggregate results
6. `6x_plot_*.R` scripts produce final figures

**Demultiplexing methods benchmarked** (all in `bin/`):
- `run_fishash.R` — the method being evaluated (0 and 10 refitting iterations)
- `run_sceptre_mixture.R` — mixture model baseline
- `run_demuxem.py` — signal thresholding (2 and 10 UMI cutoffs)
- `run_crispat_{gauss,poisgauss,poisson,negbinom}.py` — four model variants
- geomux — called via CLI (min_umi 1 and 5)
- CLEANSER — called via CLI (cs/dc modes, cutoff 0.5/0.95)

**Nextflow resource labels** in `nextflow.config`: `fishash` (2GB), `split_combine` (8GB), `data_conversion` (3GB), `sceptre`/`cleanser`/`crispat_big` (16GB base, scaled by `*MemFactor` params), `crispat_small`/`demuxem` (4–8GB base). Memory scales with `task.attempt` on retry.

**Simulation scenarios** (in `simulate/` dir, also mirrored in `bin/`):
- `varyNumGuides` — high vs low UMI/signal-to-noise regimes
- `varyMOI` — multiplicity of infection
- `varySignalNoiseCorr` — signal-noise correlation

## Key Configuration

`nextflow.config` params to tune for new runs:
- `skipCleanser`, `skipCrispatPoisson`, `skipCrispatNegBinom`, `skipSceptre` — skip slow methods
- `maxMemFactor`, `cleanserMemFactor`, `demuxemMemFactor`, `crispatBigMemFactor` — scale memory if jobs OOM
- `slurmClusterOptionsDefault` / `slurmClusterOptionsLong` — SLURM queue/QoS settings

## Known Issues

- crispat dask workers occasionally timeout; use `-resume` to retry failed tasks
- CLEANSER models must be precompiled before parallel runs (done once in `setup_env.sh`) to avoid race conditions (upstream issue: https://github.com/Gersbachlab-Bioinformatics/CLEANSER/issues/10)
- crispat uses a custom fork (`jackkamm/crispat@nohist`) with plotting disabled to avoid crashes
