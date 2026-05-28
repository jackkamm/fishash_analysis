#!/usr/bin/env Rscript

library(argparse)
library(Matrix)
library(SummarizedExperiment)

parser <- ArgumentParser(
    description = "Process individual teststats to prepare for merging them"
)
parser$add_argument(
    "--sim_label",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Identifier for the simulation replicate"
)
parser$add_argument(
    "--method",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Identifier for the assignment method"
)
parser$add_argument(
    "--sim_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--stats_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Mtx of guide assignments"
)
parser$add_argument(
    "--out_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Output Rds path"
)
args <- parser$parse_args()

sim <- readRDS(args$sim_rds)
mat_stats <- readMM(args$stats_mtx)

rownames(mat_stats) <- rownames(sim)
colnames(mat_stats) <- colnames(sim)

saveRDS(
    list(
        test_stats=mat_stats,
        sim_label=args$sim_label,
        method=args$method
    ),
    args$out_rds
)
