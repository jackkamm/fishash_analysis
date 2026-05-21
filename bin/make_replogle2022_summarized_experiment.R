#!/usr/bin/env Rscript

library(argparse)

library(magrittr)
library(Matrix)
library(SummarizedExperiment)

parser <- ArgumentParser(
    description = "Convert batch from Replogle2022 to SummarizedExperiment for benchmarking pipeline"
)

parser$add_argument(
    "--mtx_folder",
    required = TRUE,
    default = NULL,
    type = "character",
    help = "Path to folder containing mtx files"
)

parser$add_argument(
    "--sample_sheet_csv",
    required = TRUE,
    default = NULL,
    type = "character",
    help = "Path to sample sheet prepared by 22_prep_replogle2022_processing.bash"
)

parser$add_argument(
    "--idx",
    required = TRUE,
    default = NULL,
    type = "integer",
    help = "Row of sample sheet to generate summarized experiment for"
)

parser$add_argument(
    "--crispr_lib_csv",
    required = TRUE,
    default = NULL,
    type = "character",
    help="Path to k562_genomewide_library.csv.gz"
)

args <- parser$parse_args()

sample_sheet <- read.csv(args$sample_sheet_csv)
batch <- sample_sheet$sim_label[args$idx]
out_rds <- sample_sheet$path[args$idx]

prefix <- file.path(
    args$mtx_folder,
    batch
)

# Read data

mat <- readMM(paste0(prefix, "_matrix.mtx.gz"))
mat %<>% as('CsparseMatrix')

df_features <- readr::read_tsv(paste0(prefix, "_features.tsv.gz"),
                               col_names=c("ID", "Symbol", "type"))

df_barcodes <- readr::read_tsv(paste0(prefix, "_barcodes.tsv.gz"),
                               col_names="Barcode")

colnames(mat) <- df_barcodes$Barcode
rownames(mat) <- make.unique(df_features$Symbol)

mat_gex <- mat[df_features$type == "Gene Expression",]
mat_grna <- mat[df_features$type == "CRISPR Guide Capture",]

# Read crispr csv

df_k562_gwps_guides <- readr::read_csv(args$crispr_lib_csv)

df_k562_gwps_guides %<>% as.data.frame()

df_k562_gwps_guides %<>%
    dplyr::mutate(
        sgID_A_clean = gsub(",", "-", paste0(sgID_A, "_posA")),
        sgID_B_clean = gsub(",", "-", paste0(sgID_B, "_posB"))
    )

stopifnot(
    rownames(mat_grna) %in%
        with(df_k562_gwps_guides,
             c(sgID_A_clean, sgID_B_clean))
)

# Get posA and posB matrices

modmat_A <- fac2sparse(df_k562_gwps_guides[,"sgID_A_clean"])
modmat_B <- fac2sparse(df_k562_gwps_guides[,"sgID_B_clean"])

colnames(modmat_A) <- colnames(modmat_B) <- df_k562_gwps_guides$`unique sgRNA pair ID`

intersect_A <- intersect(rownames(modmat_A), rownames(mat_grna))

cntA <-t(modmat_A[intersect_A,]) %*% mat_grna[intersect_A,]

intersect_B <- intersect(rownames(modmat_B), rownames(mat_grna))

cntB <-t(modmat_B[intersect_B,]) %*% mat_grna[intersect_B,]

# Create matrix for both nonzero

stopifnot(rownames(cntA) == rownames(cntB))
stopifnot(rownames(cntA) == df_k562_gwps_guides$`unique sgRNA pair ID`)

stopifnot(colnames(cntA) == colnames(cntB))
stopifnot(colnames(cntA) == colnames(mat_grna))
stopifnot(colnames(cntA) == colnames(mat_gex))

both_nonzero <- (cntA > 0) & (cntB > 0)

# Create SummarizedExperiment

se_grna <- SummarizedExperiment(
    assays=list(
        counts=rbind(cntA, cntB),
        ground_truth=rbind(both_nonzero, both_nonzero)
    )
)

rownames(se_grna) <- with(
    df_k562_gwps_guides,
    c(sgID_A_clean, sgID_B_clean)
)

attr(se_grna, 'gex') <- SummarizedExperiment(
    assays=list(counts=mat_gex)
)

stopifnot(colnames(se_grna) == colnames(attr(se_grna, 'gex')))

# Save it

saveRDS(se_grna, out_rds)
