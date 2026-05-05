#!/usr/bin/env Rscript
library(argparse)
library(fishash)
library(Matrix)
library(SummarizedExperiment)

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
parser$add_argument(
    "--padj_cutoff",
    required = FALSE,
    type="double",
    default=0.05,
    help="FDR cutoff"
)
args <- parser$parse_args()

writeMM(
    assay(
        fishash(assay(readRDS(args$in_rds), "counts"),
                refit=args$refit,
                padj_cutoff=args$padj_cutoff,
                exclude_empty=TRUE),
        'assigned'
    ),
    args$out_mtx
)
