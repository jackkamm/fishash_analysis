source("include/env_vars.sh")

OUTS <- normalizePath(OUTS)
raw_data_dir <- file.path(OUTS, "barnyard_data/raw")
processed_data_dir <- file.path(OUTS, "barnyard_data/processed")

## Load libraries

library(magrittr)
library(ggplot2)

library(patchwork)

library(Matrix)

library(SingleCellExperiment)

library(BiocParallel)
if ("SLURM_CPUS_ON_NODE" %in% names(Sys.getenv())) {
  ncores <- as.integer(Sys.getenv()[["SLURM_CPUS_ON_NODE"]])
  BPPARAM <- MulticoreParam(ncores)
} else {
  BPPARAM <- MulticoreParam(8)
}

## Helper function to read in batches

read_liu2025cleanser_batch <- function(batch, prefix) {
    mat <- readMM(paste0(prefix, "_matrix.mtx.gz"))
    mat %<>% as('CsparseMatrix')

    df_features <- readr::read_tsv(paste0(prefix, "_features.tsv.gz"),
                                   col_names=c("ID", "Symbol", "type"),
                                   show_col_types=FALSE)

    df_barcodes <- readr::read_tsv(paste0(prefix, "_barcodes.tsv.gz"),
                                   col_names="Barcode",
                                   show_col_types=FALSE)

    df_barcodes$batch <- batch

    df_features %<>% as.data.frame()
    df_barcodes %<>% as.data.frame()

    rownames(df_barcodes) <- paste(batch, df_barcodes$Barcode, sep="_")
    rownames(df_features) <- make.unique(df_features$Symbol)

    colnames(mat) <- rownames(df_barcodes)
    rownames(mat) <- rownames(df_features)

    stopifnot(df_features$type %in% c("Gene Expression",
                                      "CRISPR Guide Capture"))

    mat_grna <- mat[df_features$type == "CRISPR Guide Capture",]
    mat_gex <- mat[df_features$type == "Gene Expression",]

    stopifnot(nrow(mat_grna) + nrow(mat_gex) == nrow(mat))

    sce <- SingleCellExperiment(
        assays=list(counts=mat_gex),
        rowData=df_features[rownames(mat_gex),],
        colData=df_barcodes,
        altExps=list(grna=SingleCellExperiment(
            assays=list(counts=mat_grna),
            rowData=df_features[rownames(mat_grna),]
        ))
    )

    stopifnot(startsWith(rowData(sce)$Symbol, 'GRCh38') |
                  startsWith(rowData(sce)$Symbol, 'mm10'))

    rowData(sce)$ref <- dplyr::if_else(
        startsWith(rowData(sce)$Symbol, 'GRCh38'),
        'homo', 'mus'
    )

    colData(sce)$homo_sum_gex <- colSums(
        counts(sce)[rowData(sce)$ref == 'homo',]
    )

    colData(sce)$mus_sum_gex <- colSums(
        counts(sce)[rowData(sce)$ref == 'mus',]
    )

    stopifnot(colData(sce)$homo_sum_gex + colData(sce)$mus_sum_gex
              == colSums(counts(sce)))

    rowData(altExp(sce))$guide_type <- dplyr::if_else(
        rownames(altExp(sce)) %in% paste0("nt_", 1:100),
        "homo_guide", "mus_guide"
    )

    colData(sce)$homo_sum_grna <- colSums(
        counts(altExp(sce))[rowData(altExp(sce))$guide_type == 'homo_guide',]
    )

    colData(sce)$mus_sum_grna <- colSums(
        counts(altExp(sce))[rowData(altExp(sce))$guide_type == 'mus_guide',]
    )

    stopifnot(colData(sce)$homo_sum_grna + colData(sce)$mus_sum_grna
              == colSums(counts(altExp(sce))))


    sce
}

## Read in the batches

prefixes <- list.files(raw_data_dir)
prefixes %<>% .[endsWith(., "_matrix.mtx.gz")]

prefixes %<>% stringr::str_match("(.*)_matrix.mtx.gz") %>% .[,2]

batch_names <- stringr::str_match(prefixes, "GSE272457_(.*)")[,2]

sce_list <- bplapply(
    1:length(prefixes),
    function(i) read_liu2025cleanser_batch(batch_names[i], file.path(raw_data_dir, prefixes[i])),
    BPPARAM=BPPARAM
)

## Variables for relabeling batches

batch_to_species_mix <- c(
    `293T_LRB100_NTlib1`="homo",
    `293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix`="mix0hr",
    `293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix`="mix72hr",
    `293T_MCH2_NTlib1`="homo",
    `293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix`="mix0hr",
    `293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix`="mix72hr",
    `NIH3T3_LRB100_NTlib2`="mus",
    `NIH3T3_MCH2_NTlib2`="mus"
)

batch_to_seqtech <- c(
    `293T_LRB100_NTlib1`="Cropseq",
    `293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix`="Cropseq",
    `293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix`="Cropseq",
    `293T_MCH2_NTlib1`="DirectCapture",
    `293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix`="DirectCapture",
    `293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix`="DirectCapture",
    `NIH3T3_LRB100_NTlib2`="Cropseq",
    `NIH3T3_MCH2_NTlib2`="DirectCapture"
)

## Helper function to merge the SCEs and give better batch names

merge_sce_list <- function(sce_list) {
    stopifnot(sapply(
        sce_list,
        function(sce) {
            all(rownames(sce) == rownames(sce_list[[1]])) &&
                all(rownames(altExp(sce)) == rownames(altExp(sce_list[[1]])))
        }
    ))

    sce_merged <- do.call(cbind, sce_list)

    # reorder batches for better plotting
    colData(sce_merged)$species_mix <- batch_to_species_mix[
        colData(sce_merged)$batch
    ]

    colData(sce_merged)$seq_tech <- batch_to_seqtech[
        colData(sce_merged)$batch
    ]

    colData(sce_merged)$batch_name <- with(
        colData(sce_merged),
        paste(species_mix, seq_tech)
    )

    sce_merged
}

## Merge the SCEs and relabel the batches

sce_merged <- merge_sce_list(sce_list)

## Split the SCEs up again (but now they have the better batch names)

metadata(sce_merged) <- list()
metadata(altExp(sce_merged)) <- list()

batches <- unique(colData(sce_merged)$batch)

sce_list <- lapply(
    batches,
    function(b) {
        sce_merged[, colData(sce_merged)$batch == b]
    }
)
names(sce_list) <- batches

## Convert the SCEs to SummarizedExperiment
## (The reason for doing this is that, while originally developing
## this script and processing the barnyard data, I used a different
## version of R, and there were problems serializing
## SingleCellExperiment between R/bioconductor versions; whereas
## serializing SummarizedExperiments works across R versions)

list_converted <- lapply(
    sce_list,
    function(x) {
        altx <- altExp(x)

        # recreate SCE from scratch to avoid unwanted KITE dependencies
        ret <- SummarizedExperiment(
            assays=list(counts=assay(altx, 'counts')),
            colData=as.data.frame(colData(altx)),
            rowData=as.data.frame(rowData(altx))
        )
        
        # HACK: Workaround so I can run the simulation nextflow
        # pipeline which requires a ground truth.
        # TODO: Refactor nextflow pipeline so methods can be run
        # without doing the confusion matrix step.
        assay(ret, 'ground_truth') <- assay(ret, 'counts') > 0

        # save the gex as well, so I can update sceptre script to use
        # it when available
        attr(ret, 'gex') <- SummarizedExperiment(
            assays=list(counts=assay(x, 'counts')),
            colData=as.data.frame(colData(x)),
            rowData=as.data.frame(
                rowData(x)[,c('ID', 'Symbol', 'type', 'ref')]
            )
        )

        ret
    }
)

## Save the individual SummarizedExperiments

dir.create(processed_data_dir)

paths <- file.path(processed_data_dir, paste0(names(list_converted), ".Rds"))

for (i in 1:length(paths)) {
    saveRDS(
        list_converted[[i]],
        paths[i]
    )
}

write.csv(
    data.frame(
        sim_label=names(list_converted),
        path=paths
    ),
    file.path(processed_data_dir, "sample_sheet.csv"),
    row.names=F, quote=F
)
