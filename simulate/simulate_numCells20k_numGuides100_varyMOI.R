library(magrittr)
library(fishash)
library(SummarizedExperiment)

out_dir <- commandArgs(trailingOnly=TRUE)[1]
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

seed <- (202512231439 %% (2^31-1))
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
        simulate_guidebender(
            n_guides=100,
            n_cells=20000,
            moi=meta_df$moi[i],
            hurdle_prob=.1,
            d_mu_cell=log(100),
            d_mu_drop=log(20),
            d_sigma_drop=.25, d_sigma_cell=.25, d_sigma_guide=.25,
            rho_alpha=5, rho_beta=45,
            return_sparse_only=TRUE
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
