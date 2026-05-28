#!/usr/bin/env Rscript

library(magrittr)
library(Matrix)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly=TRUE)

out_rds <- args[1]
sample_sheet_csv <- args[2]
in_rds_list <- args[-(1:2)]

sample_sheet <- read.csv(sample_sheet_csv)

rep_list <- lapply(
    sample_sheet$path,
    function(x) {
        list(
            simulation=readRDS(x),
            results=list()
        )
    }
)
names(rep_list) <- sample_sheet$sim_label

for (in_rds in in_rds_list) {
    curr <- readRDS(in_rds)
    rep_list[[curr$sim_label]]$results[[curr$method]] <- curr$test_stats
}

for (i in names(rep_list)) {
    metadata(rep_list[[i]]$simulation)$sim_label <- i
    rep_list[[i]]$results <- SummarizedExperiment(
        assays=rep_list[[i]]$results[sort(names(rep_list[[i]]$results))]
    )
}

saveRDS(
    rep_list,
    out_rds
)
