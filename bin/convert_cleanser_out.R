#!/usr/bin/env Rscript

library(argparse)
library(Matrix)

parser <- ArgumentParser(
    description = "Convert cleanser posterior matrix to assignments"
)

parser$add_argument(
    "--cleanser_posterior",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to cleanser posterior"
)

parser$add_argument(
    "--out_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to output matrix"
)

parser$add_argument(
    "--cutoff",
    type="double",
    default=0.5,
    help="Posterior probability cutoff"
)

args <- parser$parse_args()

# readMM doesn't work because the cleanser doesn't output the MM
# header properly. so read it in manually
mat_cleanser <- read.table(args$cleanser_posterior)

mat_cleanser_dims <- unlist(mat_cleanser[1,1:2])
names(mat_cleanser_dims) <- NULL

mat_cleanser <- mat_cleanser[-1,]

mat_cleanser <- sparseMatrix(
    i=mat_cleanser[,1], j=mat_cleanser[,2],
    x=mat_cleanser[,3],
    dims=mat_cleanser_dims,
    index1=TRUE
)

mat_cleanser <- mat_cleanser >= args$cutoff

writeMM(mat_cleanser, args$out_mtx)
