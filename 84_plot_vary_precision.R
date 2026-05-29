library(magrittr)
library(ggplot2)

library(Matrix)
library(SummarizedExperiment)

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)


## ADDITIONAL HELPER FUNCTIONS

clean_dataframe_discreteSweep <- function(df) {
    df %>%
        mutate(tuning_param=as.numeric(stringr::str_match(method, "([0-9\\.]+)$")[,2])) %>%
        mutate(
            method=case_when(
                startsWith(method, "fishash") ~ "fishash",
                startsWith(method, "cleanser_cs") ~ "cleanser_cs",
                startsWith(method, "cleanser_dc") ~ "cleanser_dc",
                startsWith(method, "geomux") ~ "geomux",
                startsWith(method, "sceptre_mixture") ~ "sceptre_mixture",
                startsWith(method, "demuxem") ~ "demuxem",
                TRUE ~ method
            )
        ) %>%
        mutate(method=factor(method, levels=method_levels)) %>%
        mutate(
            precision_nominal_bound=case_when(
                method %in% c("sceptre_mixture", "cleanser_cs", "cleanser_dc") ~ tuning_param,
                method %in% c("fishash", "geomux") ~ 1-tuning_param,
                TRUE ~ NA_real_
            )
        ) %>%
        filter(subset == "full")
}

add_precision_recall <- function(df) {
    df %>%
        mutate(
            Recall=TP/(TP+FN),
            Precision=TP/(TP+FP)
        )
}

aggregate_prec_recall <- function(df) {
    df %>%
        summarize(TN=sum(TN), FN=sum(FN), FP=sum(FP), TP=sum(TP),
                         .groups='drop') %>%
        add_precision_recall()
}

average_precision_score <- function(TN, FN, FP, TP) {
    recall <- TP / (TP + FN)
    precision <- TP / (TP + FP)

    tot <- unique(TN + FN + FP + TP)
    tot_class1 <- unique(TP + FN)

    stopifnot(length(tot) == 1)
    stopifnot(length(tot_class1) == 1)
    
    o <- order(recall)

    recall <- c(0, recall[o], 1)
    precision <- c(1,precision[o], tot_class1 / tot)

    sum(
        diff(recall) * precision[-1]
    )
}

add_avg_prec_score <- function(df) {
    df %>%
        group_by(method, sim_label, .add=TRUE) %>%
        summarize(
            AUPRC=average_precision_score(TN=TN, FN=FN, FP=FP, TP=TP),
            .groups='drop'
        )
}

discreteSweep_at_default <- function(method, tuning_param) {
    case_when(
        method=='demuxem' ~ tuning_param == 2,
        method=='sceptre_mixture' ~ tuning_param == 0.8,
        method %in% c('geomux', 'fishash') ~ tuning_param == 0.05,
        method %in% c('cleanser_cs', 'cleanser_dc') ~ tuning_param == 0.95,
        TRUE ~ as.logical(NA)
    )
}

teststats_pair2df <- function(se_pair) {
    sim <- se_pair$simulation
    res <- se_pair$results

    nonzero_mask <- (assay(sim, 'counts') > 0)
    truth_mask <- assay(sim, 'ground_truth')
    
    n_zero <- nrow(sim) * ncol(sim) - sum(nonzero_mask)
    n_zero_pos <- sum(truth_mask) - sum(truth_mask & nonzero_mask)

    ret <- list()
    for (k in names(assays(res))) {
        # we assume the teststats are 0 whenever count is 0
        stopifnot(assay(sim, 'counts')[assay(res, k) > 0] != 0)

        curr <- data.frame(
            test_stat = assay(res, k)[nonzero_mask],
            label1 = as.integer(assay(sim, 'ground_truth')[nonzero_mask])
        )
        curr$label0 <- 1 - curr$label1

        curr <- rbind(curr, data.frame(
            test_stat = 0,
            label1 = n_zero_pos,
            label0 = n_zero - n_zero_pos
        ))

        curr <- lazy_dt(curr)

        curr %<>%
            arrange(test_stat) %>%
            group_by(test_stat) %>%
            summarize(label0=sum(label0), label1=sum(label1)) %>%
            as.data.frame()

        curr$sim_label <- metadata(sim)$sim_label
        curr$method <- k

        ret[[k]] <- curr
    }

    `rownames<-`(do.call(rbind, ret), NULL)
}

add_teststat_cutoff_metrics <- function(df_teststats) {
    df_teststats %>%
        arrange(test_stat) %>%
        group_by(test_stat, .add=TRUE) %>%
        summarize(label0=sum(label0), label1=sum(label1),
                  .groups='drop_last') %>%
        arrange(desc(test_stat)) %>%
        mutate(
            TN=sum(label0) - cumsum(label0),
            FN=sum(label1) - cumsum(label1),
            FP=cumsum(label0),
            TP=cumsum(label1)
        ) %>%
        add_precision_recall()
}

# Compresses the PR curve by taking its values on a fine grid
# Speeds up plotting and reduces the size of vector graphics
compress_prc_df <- function(grouped_lazy_dt, nsteps) {
    recall_rounded <- seq(from=0, to=1, length.out=nsteps)

    grouped_lazy_dt %>%
        arrange(Recall) %>%
        # HACK use c(...) to workaround a weird dtplyr macro expansion bug
        summarize(Recall_rounded=c(recall_rounded),
                  Precision_rounded=stepfun(
                      Recall, c(Precision,0),
                      right=TRUE
                  )(recall_rounded),
                  .groups='drop') %>%
        as.data.frame() %>%
        rename(Recall=Recall_rounded, Precision=Precision_rounded)
}


## VARYING NUMBER OF GUIDES SCENARIO

prefix_high_varyPrec <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi100_snr4_endo75_varyNumGuides_varyPrecision"
)

prefix_low_varyPrec <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi20_snr1_endo25_varyNumGuides_varyPrecision"
)

prefix_high_orig <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi100_snr4_endo75_varyNumGuides"
)

prefix_low_orig <- file.path(
    OUTS,
    "results",
    "numCells20k_medUmi20_snr1_endo25_varyNumGuides"
)

# Load the results from discrete parameter sweep

df_prec_recall_varyNguides <- rbind(
    read.csv(file.path(
        prefix_high_varyPrec,
        "combined_confusion.csv"
    )) %>%
        mutate(regime="high_grna"),
    read.csv(file.path(
        prefix_low_varyPrec,
        "combined_confusion.csv"
    )) %>%
        mutate(regime="low_grna")
)

df_prec_recall_varyNguides %<>% 
    clean_dataframe_discreteSweep()

df_prec_recall_varyNguides %<>% cbind(parse_simlab_varyNumGuides(.$sim_label))

df_prec_recall_varyNguides %>%
    group_by(nguides, regime) %>%
    add_avg_prec_score() ->
    df_prec_recall_varyNguides_auprc

# Load the results from the continuous test statistics

list_teststats_highexpr <- readRDS(
    file.path(
        prefix_high_varyPrec, "combined_teststats.Rds"
    )
)

list_teststats_lowexpr <- readRDS(
    file.path(
        prefix_low_varyPrec, "combined_teststats.Rds"
    )
)

df_teststats_varyNguides <- rbind(
    lapply(
        list_teststats_highexpr,
        teststats_pair2df
    ) %>%
        rbindlist() %>%
        lazy_dt() %>%
        mutate(regime='high_grna') %>%
        as.data.frame(),
    lapply(
        list_teststats_lowexpr,
        teststats_pair2df
    ) %>%
        rbindlist() %>%
        lazy_dt() %>%
        mutate(regime='low_grna') %>%
        as.data.frame()
)

df_teststats_varyNguides %<>%
    select(sim_label) %>%
    distinct() %>%
    cbind(parse_simlab_varyNumGuides(.$sim_label)) %>%
    right_join(df_teststats_varyNguides) %>%
    mutate(method=if_else(
        method=='fishash_refit10_padj0.05',
        'fishash', method
    ))

# compute AUPRC per replicate

df_teststats_varyNguides %>%
    lazy_dt() %>%
    group_by(regime, nguides, method, sim_label) %>%
    add_teststat_cutoff_metrics() %>%
    summarize(AUPRC=average_precision_score(TN, FN, FP, TP),
              .groups='drop') %>%
    as.data.frame() ->
    df_teststats_varyNguides_auprc

# compute precision-recall for test-stats, aggregating the replicates
# within each sub-scenario
df_teststats_varyNguides %>%
    lazy_dt() %>%
    group_by(nguides, regime, method) %>%
    add_teststat_cutoff_metrics() %>%
    ungroup() %>%
    as.data.frame() ->
    df_teststats_aggPrecRecall_varyNguides

df_teststats_aggPrecRecall_varyNguides %>%
    lazy_dt() %>%
    group_by(regime, nguides, method) %>%
    compress_prc_df(1001) ->
    df_teststats_aggPrecRecall_varyNguides_compressed

# Load precision-recall from the older run to get the crispat results
# evaluated at a single point, no precision-recall tradeoff

df_prec_recall_varyNguides_crispat <- rbind(
    read.csv(file.path(
        prefix_high_orig,
        "leq20k",
        "combined_confusion.csv"
    )) %>%
        mutate(regime="high_grna"),
    read.csv(file.path(
        prefix_high_orig,
        "80k",
        "combined_confusion.csv"
    )) %>%
        mutate(regime="high_grna"),
    read.csv(file.path(
        prefix_low_orig,
        "leq20k",
        "combined_confusion.csv"
    )) %>%
        mutate(regime="low_grna"),
    read.csv(file.path(
        prefix_low_orig,
        "80k",
        "combined_confusion.csv"
    )) %>%
        mutate(regime="low_grna")
) %>%
    filter(subset=='full') %>%
    filter(startsWith(method, 'crispat'))

df_prec_recall_varyNguides_crispat %<>% cbind(parse_simlab_varyNumGuides(.$sim_label))

df_prec_recall_varyNguides_crispat %<>%
    group_by(nguides, regime, method) %>%
    summarize(TN=sum(TN), FN=sum(FN), FP=sum(FP), TP=sum(TP),
              .groups='drop') %>%
    add_precision_recall()

## Plots

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_varyPrecision_precVsRecall.pdf"),
    width=14, height=7)
bind_rows(
    df_prec_recall_varyNguides %>%
        group_by(method, regime, nguides, tuning_param, precision_nominal_bound) %>%
        summarize(TN=sum(TN), FN=sum(FN), FP=sum(FP), TP=sum(TP),
                  .groups='drop') %>%
        add_precision_recall() %>%
        mutate(curve_type="Discretized") %>%
        mutate(default=discreteSweep_at_default(method, tuning_param)),
    df_teststats_aggPrecRecall_varyNguides_compressed %>%
        mutate(curve_type="Full"),
    df_prec_recall_varyNguides_crispat %>%
        mutate(curve_type="Discretized",
               default=TRUE),
    ) %>%
    mutate(method=forcats::fct_relevel(method, method_levels)) %>%
    arrange(desc(curve_type), method) %>% {
        ggplot(., aes(x=Recall, y=Precision, color=method)) +
            geom_step(aes(lty=curve_type), direction="vh") +
            geom_point(aes(shape=default),
                       data=filter(., curve_type=='Discretized')) +
            facet_grid(regime~nguides, labeller='label_both') +
            scale_color_manual(values=method_colors) +
            scale_shape_manual(values=c(1,8)) +
            scale_linetype_manual(values=c('solid', 'dotdash')) +
            xlim(0,1) +
            ylim(0,1) +
            guides(shape=guide_legend(nrow=2), linetype=guide_legend(nrow=2),
                   color=guide_legend(nrow=2)) +
            theme_classic(base_size=16) +
            theme(axis.text.x=element_text(angle=45, hjust=1),
                  legend.box.margin = margin(t = 0, r = 0, b = 0, l = -30,
                                             unit = "pt"),
                  legend.position='bottom')
    }
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_varyPrecision_auprc_vs_f1.pdf"),
    width=14, height=7)
  set.seed(12345)
  bind_rows(
      df_teststats_varyNguides_auprc %>%
          mutate(curve_type="Full"),
      df_prec_recall_varyNguides_auprc %>%
          mutate(curve_type="Discretized")
  ) %>%
      full_join(
          df_prec_recall_varyNguides %>%
              filter(discreteSweep_at_default(method, tuning_param)) %>%
              mutate(F1=f1_score(Precision, Recall)) %>%
              select(method, sim_label, regime, nguides, F1)
      ) %>%
      .[sample(1:nrow(.)),] %>%
      ggplot(aes(x=F1, y=AUPRC, color=method, shape=curve_type)) +
      geom_point(aes(size=curve_type)) +
      scale_color_manual(values=method_colors) +
      facet_grid(regime~nguides, labeller='label_both') +
      #xlim(0,1) + ylim(0,1) +
      scale_shape_manual(values=c(1,4)) +
      scale_size_manual(values=c(2.75,1)) +
      theme_bw(base_size=16) +
      theme(legend.position='bottom')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_varyPrecision_auprc.pdf"),
    width=14, height=7)
df_prec_recall_varyNguides_auprc %>%
    ggplot(aes(x=method, y=AUPRC, color=method)) +
    geom_boxplot() +
    scale_color_manual(values=method_colors) +
    ylim(0,1) +
    facet_grid(regime~nguides) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none')
dev.off()

pdf(file.path(plot_dir, "numCells20k_varyNumGuides_varyPrecision_precVsNom.pdf"),
    width=14, height=7)
df_prec_recall_varyNguides %>%
    group_by(nguides, regime) %>%
    group_by(method, tuning_param, precision_nominal_bound, .add=TRUE) %>%
    aggregate_prec_recall() %>%
    filter(method != 'demuxem') %>%
    ggplot(aes(x=precision_nominal_bound, y=Precision, color=method)) +
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


### TODO: Update plots below!

## varyMoi scenario

df_prec_recall_varyMoi <- read.csv(
    file.path(OUTS, "results", "numCells20k_numGuides200_varyMOI_varyPrecision",
              "combined_confusion.csv")
)

df_prec_recall_varyMoi %<>% 
    clean_dataframe_discreteSweep()

df_prec_recall_varyMoi %<>% cbind(parse_simlab_varyMoi(.$sim_label))

df_prec_recall_varyMoi %>%
    group_by(moi) %>%
    group_by(method, tuning_param, precision_nominal_bound, .add=TRUE) %>%
    aggregate_prec_recall() ->
    df_prec_recall_varyMoi_sum

df_prec_recall_varyMoi %>%
    group_by(moi) %>%
    add_avg_prec_score() ->
    df_prec_recall_varyMoi_auprc

pdf(file.path(plot_dir, "numCells20k_varyMoi_varyPrecision_precVsNom.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi_sum %>%
    filter(method != 'demuxem') %>%
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

pdf(file.path(plot_dir, "numCells20k_varyMoi_varyPrecision_precVsRecall.pdf"),
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

pdf(file.path(plot_dir, "numCells20k_varyMoi_varyPrecision_auprc.pdf"),
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

pdf(file.path(plot_dir, "numCells20k_varyMoi_varyPrecision_auprc_vs_f1.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi %>%
    mutate(
        default=case_when(
            method=='demuxem' ~ tuning_param == 2,
            method=='sceptre_mixture' ~ precision_nominal_bound == 0.8,
            TRUE ~ precision_nominal_bound == 0.95
        )
    ) %>%
    filter(default) %>%
    mutate(default_f1=f1_score(Precision, Recall)) %>%
    select(method, subset, sim_label, moi, default_f1) %>%
    right_join(df_prec_recall_varyMoi_auprc) %>%
    ggplot(aes(x=default_f1, y=auprc, color=method)) +
    geom_point(shape=1) +
    scale_color_manual(values=method_colors) +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    #xlim(0,1) + ylim(0,1) +
    xlab("F1 (at default cutoffs)") +
    ylab("AUPRC (coarse-grained estimate)") +
    theme_bw(base_size=16) +
    theme(legend.position='bottom')
dev.off()
