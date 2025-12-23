library(magrittr)
library(fishash)
library(SummarizedExperiment)

out_dir <- commandArgs(trailingOnly=TRUE)[1]
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

seed <- 511082
ntimes <- 10

set.seed(seed)

meta_df <- expand.grid(
    iter=1:ntimes,
    nguides=c(20, 200, 2000, 20000, 80000)
)

meta_df %<>%
    dplyr::mutate(sim_label=sprintf("nguides_%d_iter_%d", nguides, iter)) %>%
    dplyr::mutate(path=file.path(out_dir, paste0(sim_label, ".Rds")))

for (i in 1:nrow(meta_df)) {
    saveRDS(
        simulate_guidebender2(
            n_guides=meta_df$nguides[i],
            n_cells=20000,
            moi=.3,
            hurdle_prob=.1,
            snr=4,
            count_per_cell=100,
            frac_noise_endo=.75,
            return_sparse_only=TRUE,
            chunk_cells=1000
        ),
        meta_df$path[i]
    )
}

meta_df %>%
    dplyr::select(sim_label, path) %>%
    write.csv(
        file.path(out_dir, "sample_sheet.csv"),
        row.names=F, quote=F
    )

meta_df %>%
    dplyr::filter(nguides <= 20000) %>%
    dplyr::select(sim_label, path) %>%
    write.csv(
        file.path(out_dir, "sample_sheet_leq20k.csv"),
        row.names=F, quote=F
    )

meta_df %>%
    dplyr::filter(nguides == 80000) %>%
    dplyr::select(sim_label, path) %>%
    write.csv(
        file.path(out_dir, "sample_sheet_80k.csv"),
        row.names=F, quote=F
    )
