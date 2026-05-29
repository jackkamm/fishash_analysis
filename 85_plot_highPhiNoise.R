library(magrittr)
library(ggplot2)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

## Combined plot of 20k cells with varying num guides

prefix_high <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi100_snr4_endo75_varyNumGuides_highPhiNoise"
)

prefix_low <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi20_snr1_endo25_varyNumGuides_highPhiNoise"
)

df_prec_recall <- rbind(
    read.csv(file.path(
        prefix_high, "leq20k",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="high_grna"),
    read.csv(file.path(
        prefix_high, "80k",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="high_grna"),
    read.csv(file.path(
        prefix_low, "leq20k",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="low_grna"),
    read.csv(file.path(
        prefix_low, "80k",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="low_grna")
)

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

df_prec_recall %<>% cbind(
    parse_simlab_varyNumGuides(
        .$sim_label
    )
)

df_prec_recall %<>%
    subset_methods()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_precrecall_highPhiNoise.pdf"), width=14, height=7)
set.seed(12345) # for shuffling the data point order
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    .[sample(1:nrow(.)),] %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_point(shape=8) +
    facet_grid(regime~nguides, labeller='label_both') +
    #xlim(0,1) +
    #ylim(0, 1) +
    #scale_color_manual(values=pals::okabe(3)) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_f1_highPhiNoise.pdf"), width=14, height=7)
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    ggplot(aes(x=method, y=F1, color=method)) +
    geom_boxplot() +
    facet_grid(regime~nguides, labeller='label_both') +
    ylim(0, 1) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

# top-ranked method (by median F1 score) in each scenario
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, nguides, method) %>%
    dplyr::summarize(median_f1=median(F1), .groups='drop_last') %>%
    dplyr::arrange(desc(median_f1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyNumGuides_top_per_scenario_highPhiNoise.csv"),
        row.names=F
    )

# number of simulations each method was the top one (by F1 score)
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(regime, sim_label) %>%
    dplyr::arrange(desc(F1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    .$method %>%
    table() %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyNumGuides_ntimes_top_per_sim_highPhiNoise.csv"),
        row.names=F
    )

# number of simulations each method fell below a certain precision level
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(n_below_95=sum(Precision < .95),
                     n_below_90=sum(Precision < .9)) %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyNumGuides_ntimes_precision_below_highPhiNoise.csv"),
        row.names=F
    )

## Varying MOI scenario

fname <- file.path(
    OUTS,
    "results",
    "numCells20k_numGuides200_varyMOI_highPhiNoise",
    "combined_confusion.csv"
)

df_prec_recall <- read.csv(fname)

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

df_prec_recall %<>% cbind(parse_simlab_varyMoi(.$sim_label))

df_prec_recall %<>%
    subset_methods()

pdf(file.path(plot_dir, "numCells20k_numGuides200_varyMOI_highPhiNoise_precrecall.pdf"), width=14, height=7)
set.seed(12345) # for shuffling the data point order
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    .[sample(1:nrow(.)),] %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_point(shape=8) +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    #xlim(0,1) +
    #ylim(0, 1) +
    #scale_color_manual(values=pals::okabe(3)) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_numGuides200_varyMOI_highPhiNoise_f1.pdf"), width=14, height=7)
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    ggplot(aes(x=method, y=F1, color=method)) +
    geom_boxplot() +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    ylim(0, 1) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

# top-ranked method (by median F1 score) in each scenario
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(moi, method) %>%
    dplyr::summarize(median_f1=median(F1), .groups='drop_last') %>%
    dplyr::arrange(desc(median_f1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    write.csv(
        file.path(plot_dir, "numCells20k_numGuides200_varyMOI_highPhiNoise_top_per_scenario.csv"),
        row.names=F
    )

# number of simulations each method was the top one (by F1 score)
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(sim_label) %>%
    dplyr::arrange(desc(F1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    .$method %>%
    table() %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyMOI_ntimes_top_per_sim_highPhiNoise.csv"),
        row.names=F
    )

# number of simulations each method fell below a certain precision level
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(n_below_95=sum(Precision < .95),
                     n_below_90=sum(Precision < .9)) %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyMOI_ntimes_precision_below_highPhiNoise.csv"),
        row.names=F
    )
