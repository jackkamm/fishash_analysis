#!/usr/bin/env Rscript
library(argparse)
library(SummarizedExperiment)
library(reticulate)
use_python(python=Sys.which("python"))

anndata <- import("anndata")

parser <- ArgumentParser(description = "Convert guidebender counts to h5ad")
parser$add_argument(
    "--in_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--out_h5ad",
    required = TRUE,
    default = NULL,
    type = "character",
    help="h5ad to save anndata to"
)
args <- parser$parse_args()

sim <- readRDS(args$in_rds)

adata <- anndata$AnnData(
    t(assay(sim, 'counts'))
)
adata$obs_names <- colnames(sim)
adata$var_names <- rownames(sim)

# required for some crispat models
py_run_string("r.adata.X = r.adata.X.astype(int)")
adata$obs$batch <- as.integer(0)

adata$write(args$out_h5ad)
