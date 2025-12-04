#!/usr/bin/env Rscript

library(argparse)
library(magrittr)
library(Matrix)

crispat_to_matrix <- function(path) {
    read.csv(path) %>%
        dplyr::mutate(
            row_idx=as.integer(factor(gRNA, rownames(sim))),
            col_idx=as.integer(factor(cell, colnames(sim)))
        ) %>%
        {sparseMatrix(
            i=.$row_idx,
            j=.$col_idx,
            x=rep(TRUE, nrow(.)),
            index1=TRUE,
            dims=dim(sim),
            dimnames=dimnames(sim)
        )}
}

parser <- ArgumentParser(
    description = "Convert crispat output dataframe to mtx format"
)
parser$add_argument(
    "--sim_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--crispat_out",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to crispat result"
)
parser$add_argument(
    "--out_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to save results to"
)
args <- parser$parse_args()

sim <- readRDS(args$sim_rds)
mat <- crispat_to_matrix(args$crispat_out)

writeMM(mat, args$out_mtx)
