library(magrittr)
library(ggplot2)

library(SummarizedExperiment)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

prefix <- file.path(
    OUTS,
    "results",
    "replogle2022"
)

df_prec_recall <- read.csv(file.path(
    prefix,
    "combined_confusion.csv"
))

df_prec_recall %<>%
    dplyr::filter(subset == "full") %>%
    subset_methods()

# Extract cell-level stats, used for reporting some summary stats as
# well as constructing a baseline classifier

batches <- unique(df_prec_recall$sim_label)

lapply(
    batches, function(x) {
        se <- readRDS(file.path(OUTS, "replogle2022",
                                "processed", paste0(x, ".Rds")))
        data.frame(
            sim_label=x,
            n_umi=colSums(assay(se, 'counts')),
            n_nonzero=colSums(assay(se, 'counts') > 0),
            n_true=colSums(assay(se, 'ground_truth'))
        )
    }
) %>%
    do.call(what=rbind) ->
    df_cell_stats

write.csv(df_cell_stats, file.path(OUTS, "replogle2022", "replogle2022_cell_stats.csv"))

# Add baseline of calling every guide present if it has >0 UMIs.
# (Note the proxy ground truth is a strict subset of the nonzero guides)

df_cell_stats %>%
    dplyr::group_by(sim_label) %>%
    dplyr::summarize(Precision=sum(n_true) / sum(n_nonzero)) %>%
    dplyr::mutate(Recall=1, method="baseline") %>%
    .[, c("method", "sim_label", "Precision", "Recall")] ->
    df_baseline
    
df_prec_recall %<>% .[, colnames(df_baseline)]

df_prec_recall %<>% rbind(df_baseline)

# Add F1 score

df_prec_recall %<>% dplyr::mutate(F1=f1_score(Precision, Recall))

# plots

pdf(file.path(plot_dir, "replogle2022_f1_boxplot.pdf"), width=6, height=6)
df_prec_recall %>%
    ggplot(aes(x=method, y=F1, color=method)) +
    geom_boxplot() +
    ylim(0, 1) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "replogle2022_prec_recall.pdf"), width=6, height=6)
set.seed(12345) # for shuffling the data point order
df_prec_recall %>%
    .[sample(1:nrow(.)),] %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    #geom_hline(yintercept=.95, lty='dotted') +
    #geom_point(alpha=.5, shape=1, position=position_jitterdodge(jitter.width=.1)) +
    #geom_boxplot(alpha=.5) +
    geom_point(shape=1) +
    xlim(0,1) +
    ylim(0, 1) +
    #scale_color_manual(values=pals::okabe(3)) +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

df_prec_recall %>%
    dplyr::group_by(sim_label) %>%
    dplyr::arrange(desc(F1)) %>%
    dplyr::mutate(rank=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank==1) %>%
    .$method %>%
    table() %>%
    write.csv(
        file.path(plot_dir, "replogle2022_ntimes_top_per_sim.csv"),
        row.names=F
    )

df_prec_recall %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(F1=median(F1), Precision=median(Precision), Recall=median(Recall),
                     .groups='drop') %>%
    write.csv(
        file.path(plot_dir, "replogle2022_median_metrics.csv"),
        row.names=F
    )

## Check runtime

df_trace <- read.table(
    file.path(
        OUTS,
        "results",
        "replogle2022",
        "pipeline_trace.txt"
    ),
    header=T, sep='\t'
)


matched <- stringr::str_match(df_trace$name, "(.*) \\(.*\\)$")

df_trace$method <- dplyr::if_else(
    is.na(matched[,2]),
    df_trace$name,
    matched[,2]
)

df_trace_sub <- clean_trace_df(df_trace)

df_trace_sub %<>%
    subset_methods()

pdf(file.path(plot_dir, "replogle2022_minutes_log.pdf"), width=5, height=5)
df_trace_sub %>%
    ggplot(aes(x=method, y=minutes, color=method)) +
    geom_boxplot() +
    scale_y_log10() +
    coord_flip() +
    scale_color_manual(values=method_colors) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

## save some summary stats

df_gwlib <- readr::read_csv("data/k562_genomewide_library.csv.gz")

meta_summary_stats <- data.frame(
    baseline_prec_median = median(df_baseline$Precision),
    nonzero_median = median(df_cell_stats$n_nonzero),
    n_batch = length(batches),
    n_cell = nrow(df_cell_stats),
    sgrna_pairs = nrow(df_gwlib),
    sgrna_targets = 2*nrow(df_gwlib),
    # Note: the number of grna features in the 10x *.features.tsv.gz
    # files is equal to sgrna_targets - sgrna_dups
    sgrna_dups = sum(duplicated(with(
        df_gwlib,
        c(`targeting sequence A`, `targeting sequence B`)
    )))
)

write.csv(
    meta_summary_stats,
    file.path(plot_dir, "replogle2022_meta_summary_stats.csv"),
    quote=F,
    row.names=F
)
