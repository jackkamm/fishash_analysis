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
parser$add_argument(
    "--probability_threshold",
    required = FALSE,
    type="double",
    default=0.8,
    help="Posterior probability cutoff"
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

run_sceptre_assign_grnas <- function(sceptobj) {
    tryCatch({
        assign_grnas(sceptobj,
                     method="mixture",
                     probability_threshold=args$probability_threshold,
                     # on linux, better to rely on BLAS/LAPACK instead of mclapply
                     parallel=FALSE)
    }, error = function(msg) {
        # In rare cases the initial Poisson GLM fit failed during the
        # IRLS (running into infinite weights); specifically ran into
        # this on Replogle2022 gwps batch KD8_p3_26. In that case, try
        # running with the reduced formula as suggested in
        # https://timothy-barry.github.io/sceptre-book/assign-grnas.html#sec-mixture_method
        print("Failed, trying with reduced design")
        assign_grnas(sceptobj,
                     method="mixture",
                     probability_threshold=args$probability_threshold,
                     formula_object = formula(~ log(grna_n_nonzero+1) + log(grna_n_umis+1)),
                     parallel=FALSE)
    })
}

sceptobj %<>% run_sceptre_assign_grnas()

assignments_sceptre <- get_grna_assignments(
    sceptre_object = sceptobj,
)[rownames(guide_counts),]

colnames(assignments_sceptre) <- colnames(guide_counts)

writeMM(
    assignments_sceptre,
    args$out_mtx
)
