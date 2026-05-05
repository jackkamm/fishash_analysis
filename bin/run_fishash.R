#!/usr/bin/env Rscript
library(argparse)
library(fishash)
library(Matrix)
library(SummarizedExperiment)

# Update comment below to force fishash reruns:
# 2025-11-02: currently on fishash commit: bd68f33

parser <- ArgumentParser(description = "Run fishash")
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
    help="Path to save results to"
)
parser$add_argument(
    "--refit",
    required = FALSE,
    default = 0,
    type = "integer",
    help="Number of times to refit the model"
)
args <- parser$parse_args()

writeMM(
    assay(
        fishash(assay(readRDS(args$in_rds), "counts"),
                refit=args$refit,
                exclude_empty=TRUE),
        'assigned'
    ),
    args$out_mtx
)
