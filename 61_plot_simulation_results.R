library(magrittr)
library(ggplot2)

library(SummarizedExperiment)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

## Combined plot of 20k cells with varying num guides

prefix_high <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi100_snr4_endo75_varyNumGuides"
)

prefix_low <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi20_snr1_endo25_varyNumGuides"
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

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_precrecall.pdf"), width=14, height=7)
set.seed(12345) # for shuffling the data point order
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    .[sample(1:nrow(.)),] %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_point(shape=1) +
    facet_grid(regime~nguides, labeller='label_both') +
    #xlim(0,1) +
    #ylim(0, 1) +
    #scale_color_manual(values=pals::okabe(3)) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_f1.pdf"), width=14, height=7)
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
        file.path(plot_dir, "numCells20k_varyNumGuides_top_per_scenario.csv"),
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
        file.path(plot_dir, "numCells20k_varyNumGuides_ntimes_top_per_sim.csv"),
        row.names=F
    )

# number of simulations each method fell below a certain precision level
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(n_below_95=sum(Precision < .95),
                     n_below_90=sum(Precision < .9)) %>%
    write.csv(
        file.path(plot_dir, "numCells20k_varyNumGuides_ntimes_precision_below.csv"),
        row.names=F
    )


## Varying MOI scenario

fname <- file.path(
    OUTS,
    "results",
    "numCells20k_numGuides100_varyMOI",
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

df_prec_recall %<>%
    subset_methods()

pdf(file.path(plot_dir, "numCells20k_numGuides100_varyMOI_precrecall.pdf"), width=14, height=7)
set.seed(12345) # for shuffling the data point order
df_prec_recall %>%
    dplyr::filter(subset=='full') %>%
    .[sample(1:nrow(.)),] %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_point(shape=1) +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    #xlim(0,1) +
    #ylim(0, 1) +
    #scale_color_manual(values=pals::okabe(3)) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_numGuides100_varyMOI_f1.pdf"), width=14, height=7)
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
        file.path(plot_dir, "numCells20k_numGuides100_varyMOI_top_per_scenario.csv"),
        row.names=F
    )


## Plot runtime of 20k cells low-grna vary-nguides scenario

df_trace <- rbind(
    read.table(
        file.path(
            OUTS,
            "results",
            "numCells20k_medUmi20_snr1_endo25_varyNumGuides",
            "leq20k",
            "pipeline_trace.txt"
        ),
        header=T, sep='\t'
    ) %>%
        dplyr::mutate(regime="low_grna"),
    read.table(
        file.path(
            OUTS,
            "results",
            "numCells20k_medUmi20_snr1_endo25_varyNumGuides",
            "80k",
            "pipeline_trace.txt"
        ),
        header=T, sep='\t'
    ) %>%
        dplyr::mutate(regime="low_grna"),
    read.table(
        file.path(
            OUTS,
            "results",
            "numCells20k_medUmi100_snr4_endo75_varyNumGuides",
            "leq20k",
            "pipeline_trace.txt"
        ),
        header=T, sep='\t'
    ) %>%
        dplyr::mutate(regime="high_grna"),
    read.table(
        file.path(
            OUTS,
            "results",
            "numCells20k_medUmi100_snr4_endo75_varyNumGuides",
            "80k",
            "pipeline_trace.txt"
        ),
        header=T, sep='\t'
    ) %>%
        dplyr::mutate(regime="high_grna")
)

matched <- stringr::str_match(df_trace$name, "(.*) \\(.*\\)$")

df_trace$method <- dplyr::if_else(
    is.na(matched[,2]),
    df_trace$name,
    matched[,2]
)

df_trace_sub <- clean_trace_df(df_trace)

df_trace_sub %<>% cbind(
    parse_simlab_varyNumGuides(
        .$sim_label
    )
)

df_trace_sub %<>%
    subset_methods()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_minutes.pdf"), width=10, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=minutes, color=method)) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(regime~nguides, scales='free_x', labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_minutes_log.pdf"), width=10, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=minutes, color=method)) +
    geom_boxplot() +
    scale_y_log10() +
    coord_flip() +
    facet_grid(regime~nguides, labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_subset_seconds.pdf"), width=10, height=5)
df_trace_sub %>%
    dplyr::filter(
        !startsWith(as.character(method), 'crispat'),
        !startsWith(as.character(method), 'cleanser'),
        !startsWith(as.character(method), 'sceptre')
    ) %>%
    ggplot(aes(x=method, y=seconds, color=method)) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(regime~nguides, scales='free_x', labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

df_trace_sub %>%
    dplyr::filter(method=='fishash') %>%
    .$seconds %>%
    max() %>%
    cat(file=file.path(plot_dir, "20k_fishash_max_runtime.txt"))

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_subset_peakrss.pdf"), width=10, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=peak_rss_gb, color=method)) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(regime~nguides, scales='free_x', labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_subset_peakrss_log.pdf"), width=10, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=peak_rss_gb, color=method)) +
    geom_boxplot() +
    coord_flip() +
    scale_y_log10() +
    facet_grid(regime~nguides, labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_subset_peakvmem.pdf"), width=10, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=peak_vmem_gb, color=method)) +
    geom_boxplot() +
    coord_flip() +
    facet_grid(regime~nguides, scales='free_x', labeller='label_both') +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()
