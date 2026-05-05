library(fishash)
library(SummarizedExperiment)

out_dir <- commandArgs(trailingOnly=TRUE)[1]
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

seed <- 54321
ntimes <- 3

set.seed(seed)

sim_labels <- sprintf("sim_label_%d", 1:ntimes)
paths <- file.path(out_dir, paste0(sim_labels, ".Rds"))

for (i in 1:ntimes) {
    saveRDS(
        simulate_guidebender2(
            n_guides=13,
            n_cells=21,
            moi=.3,
            hurdle_prob=.1,
            snr=1,
            count_per_cell=30,
            frac_noise_endo=.5,
            return_sparse_only=TRUE
        ),
        paths[i]
    ) 
}

write.csv(
    data.frame(
        sim_label=sim_labels,
        path=paths
    ),
    file.path(out_dir, "sample_sheet.csv"),
    row.names=F, quote=F
)
