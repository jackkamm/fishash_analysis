#!/usr/bin/env Rscript

library(argparse)
library(SummarizedExperiment)
library(magrittr)
library(Matrix)

parser <- ArgumentParser(
    description = "Convert geomux"
)

parser$add_argument(
    "--geomux_tsv",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to geomux tsv"
)

parser$add_argument(
    "--orig_rds",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to original Rds"
)

parser$add_argument(
    "--out_mtx",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to output matrix"
)

args <- parser$parse_args()


orig_cnt <- assay(readRDS(args$orig_rds), "counts")

geomux_res <- read.table(args$geomux_tsv, header=T)

geomux_res %>%
    dplyr::filter(moi > 0) %>%
    apply(1, function(x) cbind(
        unlist(stringr::str_split(x['guide_ids_original'], '\\|')),
        x['cell_id']
    )) %>%
    do.call(what=rbind) %>%
    apply(2, as.integer) %>%
    `colnames<-`(c('guide_idx0', 'cell_idx0')) %>%
    as.data.frame() %>%
    {sparseMatrix(
        i=.$guide_idx0,
        j=.$cell_idx0,
        x=rep(TRUE, nrow(.)),
        index1 = FALSE,
        dims=dim(orig_cnt),
        dimnames=dimnames(orig_cnt)
    )} ->
    geomux_assign_mat


#stopifnot(orig_cnt[geomux_assign_mat] > 0)
# same but avoids problems with logical indexing of big matrices
stopifnot(sum(geomux_assign_mat) == sum(orig_cnt*geomux_assign_mat > 0))

writeMM(geomux_assign_mat, args$out_mtx)
