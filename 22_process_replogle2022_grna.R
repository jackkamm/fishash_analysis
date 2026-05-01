library(magrittr)
library(SummarizedExperiment)

source("include/env_vars.sh")

OUTS <- normalizePath(OUTS)

replogle2022_dir <- file.path(OUTS, "replogle2022")
processed_dir <- file.path(replogle2022_dir, "split_by_batch")
dir.create(processed_dir, showWarnings=FALSE)

# The initial input is a SummarizedExperiment containing just the gRNA
# counts and some QC metrics (doesn't include the mRNA counts to save
# memory). There are separate assays for the counts from the posA
# guides and the posB guides.
#
# TODO: Put this initial input dataset on zenodo/github assets to
# ensure reproducibility
in_rds <- file.path(replogle2022_dir, "replogle2022_grna.Rds")
se_dual <- readRDS(in_rds)

# As our initial ground truth, we will say that the KO is truly
# present if each of its 2 guides has >= 1 read. The logic is that it
# is highly unlikely to sample both guides by chance in a single
# cell. Note there might be some problems with this if there are
# outlier cells with massive number of gRNA UMIs, or if a batch has a
# particular guide-pair to be widely contaminating -- however, to
# start, we will just run on everything with this imperfect ground
# truth, and then perform additional filtering downstream if needed.
assay(se_dual, 'both_nonzero') <- (
    (assay(se_dual, 'posA_counts') > 0) & (assay(se_dual, 'posB_counts') > 0)
)

# Next, unstack the posA and posB counts into a single matrix, and
# create a corresponding SummarizedExperiment
cnt_stacked <- rbind(
    assay(se_dual, 'posA_counts'),
    assay(se_dual, 'posB_counts')
)

nonzero_stacked <- rbind(
    assay(se_dual, 'both_nonzero'),
    assay(se_dual, 'both_nonzero')
)

rowdat <- rbind(
    rowData(se_dual) %>%
        as.data.frame(optional=TRUE) %>%
        dplyr::mutate(pos='posA'),
    rowData(se_dual) %>%
        as.data.frame(optional=TRUE) %>%
        dplyr::mutate(pos='posB')
)

rowdat %<>%
    dplyr::mutate(
        sgID_orig=dplyr::if_else(pos == 'posA', sgID_A, sgID_B),
        sgID_clean=dplyr::if_else(pos == 'posA', sgID_A_clean, sgID_B_clean)
    )

rowdat %<>% .[, c(
    "unique sgRNA pair ID",
    "gene", "transcript", "ensembl gene id",
    "pos", "sgID_orig", "sgID_clean"
)]

rownames(rowdat) <- rownames(cnt_stacked) <- rownames(nonzero_stacked) <- rowdat$sgID_clean

se_stacked <- SummarizedExperiment(
    assays = list(counts=cnt_stacked, ground_truth=nonzero_stacked),
    rowData = rowdat,
    colData = colData(se_dual)
)

# Finally, split the SummarizedExperiment by 10x run, and save them to rds
batches <- unique(colData(se_stacked)$`10x_Run`)

for (b in batches) {
    keep <- colData(se_stacked)$`10x_Run` == b
    se_sub <- se_stacked[,keep]
    saveRDS(se_sub, file.path(processed_dir, paste0(b, ".rds")))
}

write.csv(
    data.frame(
        sim_label=batches,
        path=file.path(processed_dir, paste0(batches, ".rds"))
    ),
    file.path(replogle2022_dir, "sample_sheet.csv"),
    row.names=F, quote=F
)
