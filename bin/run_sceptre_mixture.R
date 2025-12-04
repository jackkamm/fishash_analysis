#!/usr/bin/env Rscript

library(argparse)
library(digest)
library(magrittr)
library(Matrix)
library(splatter)
library(sceptre)
library(SummarizedExperiment)
library(SingleCellExperiment)

parser <- ArgumentParser(description = "Run sceptre mixture assignment")
parser$add_argument(
    "--in_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of SummarizedExperiment from guidebender"
)
parser$add_argument(
    "--gex_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Rds of matrix of gene expression counts"
)
parser$add_argument(
    "--out_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Mtx to save assignments to"
)
parser$add_argument(
    "--cpus",
    required = TRUE,
    default = NULL,
    type = "integer",
    help="Number of cpus for parallelization"
)
args <- parser$parse_args()

set.seed(abs(digest2int(paste0("run_sceptre_mixture.R", "_", args$in_rds))))

se <- readRDS(args$in_rds)
guide_counts <- assay(se, "counts")

gex_mat <- readRDS(args$gex_rds)

sceptobj <- import_data(
    response_matrix = gex_mat,
    grna_matrix = guide_counts,
    grna_target_data_frame = data.frame(
        grna_id = rownames(guide_counts),
        grna_target = "non-targeting"
    ),
    # I think this only affects post-processing of the guide
    # assignments; if moi low then it will annotate the doublet status
    # and assignment per cell
    moi = "high"
)

sceptobj %<>% set_analysis_parameters()

if (args$cpus == 1) {
    sceptobj %<>% assign_grnas(method="mixture",
                               parallel=FALSE)
} else {
    sceptobj %<>% assign_grnas(method="mixture",
                               parallel=TRUE,
                               n_processors=args$cpus)
}


assignments_sceptre <- get_grna_assignments(
    sceptre_object = sceptobj
)[rownames(guide_counts),]

colnames(assignments_sceptre) <- colnames(guide_counts)

writeMM(
    assignments_sceptre,
    args$out_mtx
)
