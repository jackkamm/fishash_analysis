#!/usr/bin/env Rscript

library(argparse)
library(Matrix)
library(SummarizedExperiment)

parser <- ArgumentParser(
    description = "Process individual assignment results to prepare for merging them. Computes confusion matrices for the results, and adds metadata to the assignment matrix to prepare for merging."
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
    "--assignments_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Mtx of guide assignments"
)
parser$add_argument(
    "--out_prefix",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Prefix for output files"
)
args <- parser$parse_args()

sim <- readRDS(args$sim_rds)
assn <- readMM(args$assignments_mtx)

# compute confusion matrix without coercing sparse to dense
get_confusion <- function(true, est) {
    tp <- sum(est * true)
    fn <- sum(true) - tp
    fp <- sum(est) - tp

    ## doesn't work for vector est
    #tot <- nrow(est) * ncol(est)

    # seems to work correctly even for sparseMatrix
    stopifnot(length(est) == length(true))
    tot <- length(est)
    
    tn <- tot - tp - fn - fp

    data.frame(TN=tn, FN=fn, FP=fp, TP=tp,
               Precision=tp/(tp+fp), Recall=tp/(tp+fn))
}

nz <- assay(sim, 'counts') > 0

# save confusion matrix
write.csv(
    cbind(
        data.frame(method=args$method, sim_label=args$sim_label,
                   subset=c("full", "nonzero")),
        rbind(
            get_confusion(assay(sim, 'ground_truth'), assn),
            get_confusion(assay(sim, 'ground_truth')[nz], assn[nz])
        )
    ),
    sprintf("%s_confusion.csv", args$out_prefix),
    row.names=F
)

# save the assignments matrix with dimnames and metadata
rownames(assn) <- rownames(sim)
colnames(assn) <- colnames(sim)

saveRDS(
    list(
        assignments_mat=assn,
        sim_label=args$sim_label,
        method=args$method
    ),
    sprintf("%s_matrixWithMeta.Rds", args$out_prefix)
)
