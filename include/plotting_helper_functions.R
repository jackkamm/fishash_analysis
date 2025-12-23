
method_levels <- c(
    'cleanser_cs',
    'cleanser_dc',
    'crispat_gauss',
    'crispat_poisgauss',
    'crispat_poisson',
    'crispat_negbinom',
    'sceptre_mixture',
    'demuxem',
    'geomux',
    'fishash'
)

method_colors <- pals::brewer.paired(10)
names(method_colors) <- method_levels

## Some helper functions

f1_score <- function(precision, recall) {
    2*(precision*recall) / (precision + recall)
}

get_dataframe_counts_assignments <- function(rep_list) {
    do.call(
        rbind,
        lapply(
            rep_list,
            get_counts_and_assignments
        )
    )
}

get_counts_and_assignments <- function(sim_res) {
    assay(sim_res$simulation, 'counts') %>%
        as('TsparseMatrix') %>%
        {data.frame(
            row_idx=.@i+1,
            col_idx=.@j+1,
            count=.@x
        )} %>%
        dplyr::mutate(
            cell=colnames(sim_res$simulation)[col_idx],
            grna=rownames(sim_res$simulation)[row_idx]
        ) ->
        df

    df$sim_label = metadata(sim_res$simulation)$sim_label

    idxs <- cbind(df$row_idx, df$col_idx)
    df$ground_truth <- assay(sim_res$simulation, 'ground_truth')[idxs]

    for (n in names(assays(sim_res$results))) {
        df[,n] <- assay(sim_res$results, n)[idxs]
    }

    df
}

get_dataframe_guide_abund <- function(rep_list) {
    do.call(rbind, lapply(
        rep_list,
        function(x) {
            x <- x$simulation
            cnt <- assay(x, 'counts')
            cnt_fg <- assay(x, 'counts_signal')
            cnt_bg <- cnt - cnt_fg

            cnt_fg <- rowSums(cnt_fg)
            cnt_bg <- rowSums(cnt_bg)

            df <- data.frame(
                sim_label = metadata(x)$sim_label,
                cnt_fg = cnt_fg,
                frac_fg = cnt_fg / sum(cnt_fg),
                cnt_bg = cnt_bg,
                frac_bg = cnt_bg / sum(cnt_bg)
            )
        }
    ))
}

parse_simlab_varyNumGuides <- function(simlab) {
    matched <- stringr::str_match(simlab,
                                  "nguides_(.*)_iter_(\\d+)")

    matched <- matched[,-1]
    colnames(matched) <- c("nguides", "replicate")
    matched %<>% as.data.frame()

    matched$nguides <- as.integer(matched[,'nguides'])
    matched$replicate <- as.integer(matched[,'replicate'])

    matched
}

# A few of the methods we ran with varying parameters, but we
# typically just pick one version for presentation. Namely:
#
# fishash: Ran it with 0 refitting (not shown), and 10 iterations of refitting (default, shown)
# geomux: Ran it with min umi of 1 (not shown), and min umi of 5 (default, shown)
# demuxem: Ran it with min signal 2 (shown), and min signal 10 (default, but doesn't perform well so not shown)
# cleanser: Ran it with posterior cutoff 50% (not shown), and 95% (shown)
#
# Exploratory analysis showed that the settings we picked presented
# each of the methods in the best light in most scenarios
subset_methods <- function(df) {
    df %>%
        dplyr::filter(method %in% c(
            'fishash_refit10',
            'demuxem_signal2',
            'crispat_poisson',
            'crispat_negbinom',
            'crispat_poisgauss',
            'crispat_gauss',
            'sceptre_mixture',
            'cleanser_cs_0.95',
            'cleanser_dc_0.95',
            'geomux_minumi5'
        )) %>%
        dplyr::mutate(
            method=dplyr::case_when(
                method=="fishash_refit10" ~ "fishash",
                method=="demuxem_signal2" ~ "demuxem",
                method=='geomux_minumi5' ~ 'geomux',
                method=='cleanser_cs_0.95' ~ 'cleanser_cs',
                method=='cleanser_dc_0.95' ~ 'cleanser_dc',
                TRUE ~ method
            )
        ) %>%
        dplyr::mutate(
            method=factor(method, method_levels)
        )
}

# for plotting f1, recall, precision on subset of methods
get_long_subset_df <- function(df) {
    df %>%
        dplyr::filter(subset=='full') %>%
        #dplyr::select(method, nguides, Precision, Recall, F1) %>%
        tidyr::pivot_longer(c(Precision, Recall, F1),
                            names_to='stat', values_to='value') %>%
        subset_methods()
}

mem_string_to_gb <- function(mem_str) {
    unit <- c(GB=1, MB=1/1024)
    matched <- stringr::str_match(mem_str, "(.*) (.*)")
    as.numeric(matched[,2]) * unit[matched[,3]]
}

# for extracting trace report and plot timings
clean_trace_df <- function(df_trace) {
    df_trace %>%
        dplyr::filter(status %in% c('COMPLETED', 'CACHED')) %>%
        dplyr::filter(startsWith(method, 'run_')) %>%
        dplyr::mutate(sim_label=stringr::str_match(tag, '.*simlab=(.*)')[,2]) %>%
        dplyr::mutate(method=stringr::str_match(method, 'run_(.*)')[,2]) %>%
        dplyr::mutate(method=dplyr::if_else(
            method=='cleanser',
            paste0('cleanser_', stringr::str_match(tag, 'seqtech=(.*),simlab=.*')[,2], "_0.95"),
            method
        )) %>%
        dplyr::mutate(method=dplyr::if_else(
            method=='fishash',
            paste0('fishash_refit', stringr::str_match(tag, 'refit=(.*),simlab=.*')[,2]),
            method
        )) %>%
        dplyr::mutate(method=dplyr::if_else(
            method=='geomux',
            paste0('geomux_minumi', stringr::str_match(tag, 'minumi=(.*),simlab=.*')[,2]),
            method
        )) %>%
        dplyr::mutate(method=dplyr::if_else(
            method=='demuxem',
            paste0('demuxem_signal', stringr::str_match(tag, 'minsignal=(.*),simlab=.*')[,2]),
            method
        )) %>%
        dplyr::mutate(seconds=as.numeric(stringr::str_match(
            lubridate::seconds(lubridate::as.period(toupper(realtime))),
            "(.*)S"
        )[,2])) %>%
        dplyr::mutate(minutes=seconds / 60) %>%
        dplyr::mutate(
            peak_rss_gb = mem_string_to_gb(peak_rss),
            peak_vmem_gb = mem_string_to_gb(peak_vmem)
        )
}
