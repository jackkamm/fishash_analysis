#!/usr/bin/env Rscript

library(argparse)
library(SummarizedExperiment)
library(Matrix)

parser <- ArgumentParser(
    description = "Convert guidebender counts to mtx"
)
parser$add_argument(
    "--in_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--out_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="mtx to save counts to"
)
args <- parser$parse_args()

sim <- readRDS(args$in_rds)

writeMM(assay(sim, 'counts'), args$out_mtx)
