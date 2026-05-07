library(magrittr)
library(ggplot2)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

df_prec_recall_varyMoi <- read.csv(
    file.path(OUTS, "results", "numCells20k_numGuides200_varyMOI_varyPrecision",
              "combined_confusion.csv")
)

df_prec_recall_varyMoi %<>%
    dplyr::mutate(
        p=as.numeric(stringr::str_match(method, ".*(0.\\d+)$")[,2])
    ) %>%
    dplyr::mutate(
        method=dplyr::case_when(
            startsWith(method, "fishash") ~ "fishash",
            startsWith(method, "cleanser_cs") ~ "cleanser_cs",
            startsWith(method, "cleanser_dc") ~ "cleanser_dc",
            startsWith(method, "geomux") ~ "geomux",
            startsWith(method, "sceptre_mixture") ~ "sceptre_mixture",
            TRUE ~ method
        )
    ) %>%
    dplyr::mutate(method=factor(method, levels=method_levels))

df_prec_recall_varyMoi %<>% cbind(parse_simlab_varyMoi(.$sim_label))

df_prec_recall_varyMoi %>%
    dplyr::group_by(method, subset, p, moi) %>%
    dplyr::summarize(TN=sum(TN), FN=sum(FN), FP=sum(FP), TP=sum(TP),
                     .groups='drop') %>%
    dplyr::mutate(
        Precision=TP/(TP+FP), Recall=TP/(TP+FN)
    ) ->
    df_prec_recall_varyMoi_sum

pdf(file.path(plot_dir, "numCells20k_numGuides200_varyMOI_varyPrecision.pdf"),
    width=14, height=7)
df_prec_recall_varyMoi_sum %>%
    dplyr::filter(subset=="full") %>%
    #dplyr::filter(subset=="nonzero") %>%
    ggplot(aes(x=Recall, y=Precision, color=method)) +
    geom_line() +
    geom_point() +
    facet_wrap(~moi, labeller='label_both', nrow=2) +
    scale_color_manual(values=method_colors) +
    #xlim(0,1) +
    #ylim(0,1) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
dev.off()
