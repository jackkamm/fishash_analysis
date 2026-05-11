library(magrittr)
library(ggplot2)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

## varyNguides scenario

df_prec_recall_varyNguides <- rbind(
    read.csv(file.path(
        OUTS,
        "results",
        "numCells20k_medUmi100_snr4_endo75_varyNumGuides_varyPrecision",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="high_grna"),
    read.csv(file.path(
        OUTS,
        "results",
        "numCells20k_medUmi20_snr1_endo25_varyNumGuides_varyPrecision",
        "combined_confusion.csv"
    )) %>%
        dplyr::mutate(regime="low_grna")
)

df_prec_recall_varyNguides %<>% 
    clean_dataframe_varyPrec()

df_prec_recall_varyNguides %<>% cbind(parse_simlab_varyNumGuides(.$sim_label))

df_prec_recall_varyNguides %>%
    dplyr::group_by(nguides, regime) %>%
    aggregate_prec_recall() ->
    df_prec_recall_varyNguides_sum

df_prec_recall_varyNguides %>%
    dplyr::group_by(nguides, regime) %>%
    add_avg_prec_score() ->
    df_prec_recall_varyNguides_auprc

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_precVsNom.pdf"),
    width=14, height=7)
df_prec_recall_varyNguides_sum %>%
    dplyr::filter(method != 'demuxem') %>%
    ggplot(aes(x=precision_nominal_bound, y=Precision_agg, color=method)) +
    geom_line() +
    geom_point() +
    geom_abline(lty='dotted') +
    facet_grid(regime~nguides, labeller='label_both') +
    scale_color_manual(values=method_colors) +
    xlim(0,1) +
    ylim(0,1) +
    xlab("Precision (nominal lower bound)") +
    ylab("Precision (actual)") +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()


pdf(file.path(plot_dir, "numCells20k_varyNumGuides_precVsRecall.pdf"),
    width=14, height=7)
df_prec_recall_varyNguides_sum %>%
    ggplot(aes(y=Precision_agg, x=Recall_agg, color=method)) +
    geom_step(direction="vh") +
    geom_point() +
    facet_grid(regime~nguides, labeller='label_both') +
    scale_color_manual(values=method_colors) +
    xlim(0,1) +
    ylim(0,1) +
    xlab("Recall") +
    ylab("Precision") +
    coord_flip() +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_auprc.pdf"),
    width=14, height=7)
df_prec_recall_varyNguides_auprc %>%
    ggplot(aes(x=method, y=auprc, color=method)) +
    geom_boxplot() +
    scale_color_manual(values=method_colors) +
    ylim(0,1) +
    facet_grid(regime~nguides) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

## varyMoi scenario

df_prec_recall_varyMoi <- read.csv(
    file.path(OUTS, "results", "numCells20k_numGuides200_varyMOI_varyPrecision",
              "combined_confusion.csv")
)

df_prec_recall_varyMoi %<>% 
    clean_dataframe_varyPrec()

df_prec_recall_varyMoi %<>% cbind(parse_simlab_varyMoi(.$sim_label))

df_prec_recall_varyMoi %>%
    dplyr::group_by(moi) %>%
    aggregate_prec_recall() ->
    df_prec_recall_varyMoi_sum

df_prec_recall_varyMoi %>%
    dplyr::group_by(moi) %>%
    add_avg_prec_score() ->
    df_prec_recall_varyMoi_auprc

pdf(file.path(plot_dir, "numCells20k_varyMoi_precVsNom.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi_sum %>%
    dplyr::filter(method != 'demuxem') %>%
    ggplot(aes(x=precision_nominal_bound, y=Precision_agg, color=method)) +
    geom_line() +
    geom_point() +
    geom_abline(lty='dotted') +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    scale_color_manual(values=method_colors) +
    xlim(0,1) +
    ylim(0,1) +
    xlab("Precision (nominal lower bound)") +
    ylab("Precision (actual)") +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyMoi_precVsRecall.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi_sum %>%
    ggplot(aes(y=Precision_agg, x=Recall_agg, color=method)) +
    geom_step(direction="vh") +
    geom_point() +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    scale_color_manual(values=method_colors) +
    xlim(0,1) +
    ylim(0,1) +
    xlab("Recall") +
    ylab("Precision") +
    coord_flip() +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyMoi_auprc.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi_auprc %>%
    ggplot(aes(x=method, y=auprc, color=method)) +
    geom_boxplot() +
    scale_color_manual(values=method_colors) +
    ylim(0,1) +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
