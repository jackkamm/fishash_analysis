library(magrittr)
library(ggplot2)

library(SummarizedExperiment)

source("include/env_vars.sh")
source("include/plotting_helper_functions.R")

plot_dir <- file.path(OUTS, "plots")
dir.create(plot_dir)

## Load the full results


rep_list <- readRDS(
    file.path(
        OUTS,
        "results",
        "barnyard_data",
        "combined_results.Rds"
    )
)

## Get SE of GEX


se_gex <- do.call(
    cbind,
    lapply(
        rep_list,
        function(x) attr(x$simulation, 'gex')
    )
)

## Add some QC stats to coldata


cnt <- assay(se_gex, 'counts')

colData(se_gex)$sum_gex <- with(
    as.data.frame(colData(se_gex)),
    homo_sum_gex + mus_sum_gex
)

colData(se_gex)$ratio_gex <- with(
    as.data.frame(colData(se_gex)),
    (homo_sum_gex+1) / (mus_sum_gex+1)
)

colData(se_gex)$frac_homo <- with(
    as.data.frame(colData(se_gex)),
    (homo_sum_gex) / (homo_sum_gex + mus_sum_gex)
)

colData(se_gex)$homo_features_gex <- colSums(
    cnt[rowData(se_gex)$ref == 'homo',] > 0
)

colData(se_gex)$mus_features_gex <- colSums(
    cnt[rowData(se_gex)$ref == 'mus',] > 0
)

colData(se_gex)$features_gex <- colSums(cnt > 0)

stopifnot(with(as.data.frame(colData(se_gex)),
               features_gex == homo_features_gex + mus_features_gex))

rowData(se_gex)$is_mito <- (
    startsWith(rowData(se_gex)$Symbol, "GRCh38_MT-") |
        startsWith(rowData(se_gex)$Symbol, "mm10___mt-")
)

colData(se_gex)$homo_mito_sum <- colSums(
    cnt[with(rowData(se_gex), is_mito & ref == 'homo'),]
)

colData(se_gex)$mus_mito_sum <- colSums(
    cnt[with(rowData(se_gex), is_mito & ref == 'mus'),]
)

colData(se_gex)$mito_sum <- colSums(
    cnt[rowData(se_gex)$is_mito,]
)

stopifnot(with(colData(se_gex),
               mito_sum == homo_mito_sum + mus_mito_sum))

colData(se_gex)$mito_frac <- with(colData(se_gex),
                                      mito_sum / sum_gex)

# same cutoff for empty droplets as in CLEANSER Fig 2I
colData(se_gex)$is_empty <- colData(se_gex)$sum_gex <= 100

# from original cleanser paper
colData(se_gex)$orig_qc_pass <- with(
    colData(se_gex),
    (mito_frac < .15) &
        (features_gex <= 6000) &
        (sum_gex <= 20000) &
        (1500 <= features_gex) &
        (3500 <= sum_gex) &
        (frac_homo < .1 | frac_homo > .9)
)

# cells already filtered as high expression, but which I am keeping
# for doublet analysis
colData(se_gex)$high_expression <- with(
    colData(se_gex),
    (features_gex > 6000) | (sum_gex > 20000)
)

colData(se_gex)$has_homo <- with(
    colData(se_gex),
    (frac_homo >= .1) & !is_empty
)

colData(se_gex)$has_mus <- with(
    colData(se_gex),
    (frac_homo <= .9) & !is_empty
)

colData(se_gex)$species <- with(
    colData(se_gex),
    dplyr::case_when(
        has_homo & has_mus ~ 'both',
        has_homo ~ 'homo',
        has_mus ~ 'mus',
        TRUE ~ 'neither'
    )
)

## Get matrices of number of assigned human and mouse guides per cell


guide_type <- rowData(rep_list[[1]]$simulation)$guide_type
stopifnot(sapply(
    rep_list, function(x) rowData(x$simulation)$guide_type == guide_type
))

get_assigned_grna_cnt <- function(gtype) {
    do.call(
        rbind,
        lapply(
            rep_list,
            function(x) {
                sapply(assays(x$results),
                       function(mat) {
                           colSums(mat[guide_type == gtype,])
                       })
            }
        )
    )
}

n_assigned_mus <- t(get_assigned_grna_cnt('mus_guide'))
n_assigned_homo <- t(get_assigned_grna_cnt('homo_guide'))

## Create the SummarizedExperiment


se_assignment <- SummarizedExperiment(
    assays = list(n_assigned_homo = n_assigned_homo,
                  n_assigned_mus = n_assigned_mus),
    colData = colData(se_gex)
)

assay(se_assignment, 'has_homo') <- assay(se_assignment, 'n_assigned_homo') > 0
assay(se_assignment, 'has_mus') <- assay(se_assignment, 'n_assigned_mus') > 0

assigned_type <- matrix('neither', nrow=nrow(se_assignment), ncol=ncol(se_assignment))
assigned_type[assay(se_assignment, 'has_mus')] <- 'mus'
assigned_type[assay(se_assignment, 'has_homo')] <- 'homo'
assigned_type[assay(se_assignment, 'has_homo') &
                  assay(se_assignment, 'has_mus')] <- 'both'

assay(se_assignment, 'assigned_type', withDimnames=FALSE) <- assigned_type

assigned_type_binarized <- assay(se_assignment, 'has_homo')
assigned_type_binarized[assigned_type == 'both'] <- NA
assigned_type_binarized[assigned_type == 'neither'] <- NA

assay(se_assignment, 'assigned_type_binarized', withDimnames=FALSE) <- assigned_type_binarized

batches <- unique(colData(se_assignment)$batch)

list_se_assignment <- lapply(
    batches,
    function(b) {
        se_assignment[, colData(se_assignment)$batch == b]
    }
)
names(list_se_assignment) <- batches

## Function to compute binary accuracy metrics


binary_accuracy_metrics <- function(pred, truth) {
    tp <- sum(pred & truth, na.rm=T)
    fp <- sum(pred & !truth, na.rm=T)
    tn <- sum(!pred & !truth, na.rm=T)
    fn <- sum(!pred & truth, na.rm=T)

    naTrue <- sum(is.na(pred) & truth)
    naFalse <- sum(is.na(pred) & !truth)

    tot <- tp + fp + tn + fn + naTrue + naFalse
    tot_pos <- tp + fp
    tot_neg <- tn + fn
    tot_true <- tp + fn + naTrue
    tot_false <- tn + fp + naFalse

    stopifnot(tot == length(truth))

    acc_true <- tp / tot_true
    acc_false <- tn / tot_false

    acc_stderr <- tot_true * acc_true * (1-acc_true) + tot_false * acc_false * (1 - acc_false)
    acc_stderr <- acc_stderr / tot^2
    acc_stderr <- sqrt(acc_stderr)

    c(
        TP=tp,
        FP=fp,
        TN=tn,
        FN=fn,
        Precision=tp / tot_pos,
        Recall=tp / tot_true,
        Specificity = tn / tot_false,
        Accuracy=(tp + tn) / tot,
        Accuracy_stderr=acc_stderr,
        F1=2*tp / (tot_pos + tot_true),
        Antiprecision=tn / tot_neg,
        AntiF1 = 2*tn / (tot_neg + tot_false)
    )
}

## Compute binary accuracy metrics on each batch and original QC


bin_metrics_origqc <- lapply(
    list_se_assignment,
    function(x) {
        x <- x[, colData(x)$orig_qc_pass]
        ret <- apply(assay(x, 'assigned_type_binarized'), 1,
              function(row) {
                  binary_accuracy_metrics(row, colData(x)$has_homo)
              })
        ret <- t(ret)
        ret <- as.data.frame(ret)
        ret %<>% tibble::rownames_to_column('method')
        ret <- cbind(
            dplyr::distinct(as.data.frame(
                colData(x)[, c('batch', 'species_mix', 'seq_tech', 'batch_name')]
            )),
            ret
        )
        ret
    }
)

bin_metrics_origqc %<>% do.call(what=rbind)

## Plot accuracy

pdf(file.path(plot_dir, "barnyard_accuracy_barplot.pdf"))
bin_metrics_origqc %>%
    subset_methods() %>%
    dplyr::filter(!(species_mix %in% c('homo', 'mus'))) %>%
    ggplot(aes(x=method, y=Accuracy, fill=method)) +
    geom_col() +
    scale_fill_manual(values=method_colors) +
    facet_grid(seq_tech~species_mix) +
    theme_bw(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position='none')
dev.off()


## xtable with accuracy and stderr


library(xtable)

bin_metrics_origqc %>%
    subset_methods() %>%
    dplyr::group_by(batch_name) %>%
    dplyr::mutate(is_best=Accuracy==max(Accuracy)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!(species_mix %in% c('homo', 'mus'))) %>%
    dplyr::mutate(Accuracy_str=sprintf("$%0.4f \\pm %.04f$", Accuracy, Accuracy_stderr)) %>%
    dplyr::mutate(Accuracy_str=dplyr::if_else(
        is_best,
        paste0("\\boldmath", Accuracy_str),
        Accuracy_str
    )) %>%
    dplyr::mutate(batch_name=paste(seq_tech, species_mix)) %>%
    dplyr::select(method, Accuracy_str, batch_name) %>%
    `rownames<-`(NULL) %>%
    tidyr::pivot_wider(names_from=method, values_from=Accuracy_str) %>%
    dplyr::arrange(batch_name) %>%
    tibble::column_to_rownames("batch_name") %>%
    as.matrix() %>%
    t() %>%
    xtable() ->
    xtab

print(xtab, file=file.path(plot_dir, "barnyard_accuracy_table.tex"),
      floating=FALSE,
      sanitize.text.function=as.is,
      sanitize.rownames.function=function(x) sanitize(x, type="latex"))

## Plot num assigned of each type per method

pdf(file.path(plot_dir, "barnyard_assignedtype_barplot.pdf"), width=9.5, height=6)
cbind(
    colData(se_assignment)[, c('seq_tech', 'species_mix', 'orig_qc_pass', 'species')],
    t(assay(se_assignment, 'assigned_type'))
) %>%
    as.data.frame() %>%
    dplyr::filter(orig_qc_pass, !(species_mix %in% c('homo', 'mus'))) %>%
    tidyr::pivot_longer(
        dplyr::all_of(rownames(se_assignment)),
        names_to='method', values_to='assigned_type'
    ) %>%
    subset_methods() %>%
    ggplot(aes(x=method, fill=assigned_type)) +
    geom_bar() +
    facet_grid(species ~ seq_tech + species_mix) +
    theme_bw(base_size=16) +
    scale_fill_manual(values=pals::okabe(4)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(legend.position='bottom')
dev.off()
