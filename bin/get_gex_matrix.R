#!/usr/bin/env Rscript

library(argparse)
library(digest)
library(magrittr)
library(Matrix)
library(splatter)
library(sceptre)
library(SummarizedExperiment)
library(SingleCellExperiment)

parser <- ArgumentParser(description = "Create or extract gene expression counts from guidebender SummarizedExperiment")
parser$add_argument(
    "--in_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--out_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds path to save sceptre object"
)
args <- parser$parse_args()

se <- readRDS(args$in_rds)

guide_counts <- assay(se, "counts")

# get or simulate the GEX matrix
if (is.null(attr(se, 'gex'))) {
    set.seed(abs(digest2int(paste0("create_sceptre_obj.R", "_", args$in_rds))))

    # FIXME: It would be better to simulate this upfront with the initial
    # simulation, so that GEX and gRNA could be better correlated and
    # sceptre could make use of that information in its mixture model
    sce <- splatSimulate(
        batchCells = ncol(guide_counts),
        )
    colnames(sce) <- colnames(guide_counts)
    gex_mat <- counts(sce)
} else {
    gex_mat <- assay(attr(se, 'gex'), 'counts')
}

saveRDS(gex_mat, args$out_rds)
