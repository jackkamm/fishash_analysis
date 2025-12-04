#!/usr/bin/env Rscript

library(magrittr)
library(ggplot2)
library(Matrix)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly=TRUE)

out_prefix <- args[1]
in_csv_list <- args[-1]

combined_confusion <- do.call(
    rbind,
    lapply(in_csv_list, read.csv)
)

# full should always be bigger than nonzero
combined_confusion %>%
    dplyr::select(sim_label, method, subset, TP, FP, TN, FN) %>%
    tidyr::pivot_longer(c(TP, FP, TN, FN),
                        names_to="type", values_to="Freq") %>%
    tidyr::pivot_wider(names_from=subset, values_from=Freq) %>%
    {stopifnot(.$full >= .$nonzero)}

# no reasonable method should ever assign a zero
combined_confusion %>%
    dplyr::select(sim_label, method, subset, TP, FP) %>%
    tidyr::pivot_longer(c(TP, FP),
                        names_to="type", values_to="Freq") %>%
    tidyr::pivot_wider(names_from=subset, values_from=Freq) %>%
    {stopifnot(.$full == .$nonzero)}

write.csv(combined_confusion, sprintf("%s_confusion.csv", out_prefix),
          row.names=F)

png(paste0(out_prefix, '_precision.png'))
combined_confusion %>%
    ggplot(aes(x=method, y=Precision)) +
    geom_boxplot(alpha=.5) +
    geom_point(shape=1, alpha=.5, position=position_jitter(width=.25)) +
    facet_wrap(~subset, nrow=2) +
    ylim(0,1) +
    theme_classic(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

png(paste0(out_prefix, '_recall.png'))
combined_confusion %>%
    ggplot(aes(x=method, y=Recall)) +
    geom_boxplot(alpha=.5) +
    geom_point(shape=1, alpha=.5, position=position_jitter(width=.25)) +
    facet_wrap(~subset, nrow=2) +
    ylim(0,1) +
    theme_classic(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
