library(magrittr)
library(fishash)
library(SummarizedExperiment)

out_dir <- commandArgs(trailingOnly=TRUE)[1]
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

# seed with ISO8601 datetime at script creation
seed <- (202601050822 %% (2^31 - 1))
ntimes <- 20

set.seed(seed)

meta_df <- data.frame(
    endo_exo_corr=c("high", "mid", "low"),
    endo_shape_uniform=c(0, 0, 1),
    endo_shape_sum=c(1e6,1,1)
) %>%
    dplyr::cross_join(data.frame(
        iter=1:ntimes
    )) %>%
    dplyr::mutate(sim_label=sprintf("corr_%s_unif_%f_sum_%g_iter_%d", endo_exo_corr, endo_shape_uniform, endo_shape_sum, iter)) %>%
    dplyr::mutate(path=file.path(out_dir, paste0(sim_label, ".Rds")))

for (i in 1:nrow(meta_df)) {
    saveRDS(
        simulate_guidebender2(
            n_guides=20,
            n_cells=20000,
            moi=.3,
            hurdle_prob=.1,
            snr=4,
            count_per_cell=100,
            frac_noise_endo=1,
            endo_shape_flat=meta_df$endo_shape_uniform[i],
            endo_shape_sum=meta_df$endo_shape_sum[i],
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
