library(magrittr)
library(ggplot2)

library(SummarizedExperiment)

out_dir <- "outs"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir)

source("40_plotting_helper_functions.R")

## Combined plot of varying guides scenarios, 2000 cells

fname_high <- file.path(
    out_dir,
    "results",
    "numCells2000_medUmi100_snr4_endo75_varyNumGuides",
    "combined_confusion.csv"
)

fname_low <- file.path(
    out_dir,
    "results",
    "numCells2000_medUmi20_snr1_endo25_varyNumGuides",
    "combined_confusion.csv"
)

df_prec_recall <- rbind(
    read.csv(fname_high) %>%
        dplyr::mutate(regime="high_grna"),
    read.csv(fname_low) %>%
        dplyr::mutate(regime="low_grna")
)

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

df_prec_recall %<>% cbind(
    parse_simlab_varyNumGuides(
        .$sim_label
    )
)

pdf(file.path(plot_dir, "2k_varyguides_stats.pdf"), width=14, height=7)
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::select(method, regime, nguides, Precision, Recall, F1) %>%
    tidyr::pivot_longer(c(Precision, Recall, F1),
                        names_to='stat', values_to='value') %>%
    ggplot(aes(x=method, y=value, color=stat)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_boxplot() +
    facet_grid(regime~nguides, labeller='label_both') +
    ylim(0, 1) +
    scale_color_manual(values=pals::okabe(3)) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

# top-ranked method (by median F1 score) in each scenario
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, nguides, method) %>%
    dplyr::summarize(median_f1=median(F1), .groups='drop_last') %>%
    dplyr::arrange(desc(median_f1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    write.csv(
        file.path(plot_dir, "2k_varyguides_top_per_scenario.csv")
    )

# number of simulations each method was the top one (by F1 score)
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, sim_label) %>%
    dplyr::arrange(desc(F1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    .$method %>%
    table() %>%
    write.csv(
        file.path(plot_dir, "2k_varyguides_ntimes_top_per_sim.csv")
    )

# number of simulations each method fell below a certain precision level
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(n_below_95=sum(Precision < .95),
                     n_below_90=sum(Precision < .9)) %>%
    write.csv(
        file.path(plot_dir, "2k_varyguides_ntimes_precision_below.csv")
    )


## Combined plot of 20k cells with varying num guides

fname_high <- file.path(
    out_dir,
    "results",
    "numCells20k_medUmi100_snr4_endo75_varyNumGuides",
    "combined_confusion.csv"
)

fname_low <- file.path(
    out_dir,
    "results",
    "numCells20k_medUmi20_snr1_endo25_varyNumGuides",
    "combined_confusion.csv"
)

df_prec_recall <- rbind(
    read.csv(fname_high) %>%
        dplyr::mutate(regime="high_grna"),
    read.csv(fname_low) %>%
        dplyr::mutate(regime="low_grna")
)

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

df_prec_recall %<>% cbind(
    parse_simlab_varyNumGuides(
        .$sim_label
    )
)

pdf(file.path(plot_dir, "20k_varyguides_stats.pdf"), width=14, height=7)
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::select(method, regime, nguides, Precision, Recall, F1) %>%
    tidyr::pivot_longer(c(Precision, Recall, F1),
                        names_to='stat', values_to='value') %>%
    ggplot(aes(x=method, y=value, color=stat)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_boxplot() +
    facet_grid(regime~nguides, labeller='label_both') +
    ylim(0, 1) +
    scale_color_manual(values=pals::okabe(3)) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

# top-ranked method (by median F1 score) in each scenario
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, nguides, method) %>%
    dplyr::summarize(median_f1=median(F1), .groups='drop_last') %>%
    dplyr::arrange(desc(median_f1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    write.csv(
        file.path(plot_dir, "20k_varyguides_top_per_scenario.csv")
    )

# number of simulations each method was the top one (by F1 score)
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, sim_label) %>%
    dplyr::arrange(desc(F1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    .$method %>%
    table() %>%
    write.csv(
        file.path(plot_dir, "20k_varyguides_ntimes_top_per_sim.csv")
    )

# number of simulations each method fell below a certain precision level
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(n_below_95=sum(Precision < .95),
                     n_below_90=sum(Precision < .9)) %>%
    write.csv(
        file.path(plot_dir, "20k_varyguides_ntimes_precision_below.csv")
    )


## Varying MOI scenario

fname <- file.path(
    out_dir,
    "results",
    "numCells2000_numGuides100_varyMOI",
    "combined_confusion.csv"
)

df_prec_recall <- read.csv(fname)

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

matched <- stringr::str_match(df_prec_recall$sim_label,
                              "moi_(.*)_iter_(\\d+)")

matched <- matched[,-1]
colnames(matched) <- c("moi", "replicate")

df_prec_recall$moi <- as.numeric(matched[,'moi'])
df_prec_recall$replicate <- as.integer(matched[,'replicate'])

pdf(file.path(plot_dir, "2k_varymoi_stats.pdf"), width=14, height=7)
df_prec_recall %>%
    get_long_subset_df() %>%
    ggplot(aes(x=method, y=value, color=stat)) +
    #geom_point(alpha=.5, shape=1, position=position_jitter(width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_boxplot() +
    scale_color_manual(values=pals::okabe(3)) +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

# top-ranked method (by median F1 score) in each scenario
df_prec_recall %>%
    subset_methods() %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(moi, method) %>%
    dplyr::summarize(median_f1=median(F1), .groups='drop_last') %>%
    dplyr::arrange(desc(median_f1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    write.csv(
        file.path(plot_dir, "2k_varymoi_top_per_scenario.csv")
    )
