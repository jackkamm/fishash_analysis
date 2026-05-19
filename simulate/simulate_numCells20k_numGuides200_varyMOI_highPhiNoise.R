library(magrittr)
library(fishash)
library(SummarizedExperiment)

out_dir <- commandArgs(trailingOnly=TRUE)[1]
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

# seed with ISO8601 datetime at script creation
seed <- (202605181136 %% (2^31-1))
ntimes <- 10

set.seed(seed)

meta_df <- expand.grid(
    iter=1:ntimes,
    moi=c(.1, .3, .5, 1, 2, 3, 5, 10)
)

meta_df %<>%
    dplyr::mutate(sim_label=sprintf("moi_%f_iter_%d", moi, iter)) %>%
    dplyr::mutate(path=file.path(out_dir, paste0(sim_label, ".Rds")))

for (i in 1:nrow(meta_df)) {
    saveRDS(
        simulate_guidebender2(
            n_guides=200,
            n_cells=20000,
            moi=meta_df$moi[i],
            hurdle_prob=.1,
            snr=4,
            count_per_cell=100,
            frac_noise_endo=.75,
            return_sparse_only=TRUE,
            chunk_cells=1000,
            Phi_cell = 1,
            Phi_noise = 1
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
