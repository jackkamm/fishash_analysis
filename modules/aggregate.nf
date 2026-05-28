// Computes confusion matrices and adds metadata to assignment matrices
process process_assignments {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method_id), path(sim_rds), path(result_mtx)

    output:
    path "${sim_label}__${method_id}_confusion.csv", emit: confusion
    path "${sim_label}__${method_id}_matrixWithMeta.Rds", emit: matrix

    script:
    """
    process_assignments.R \
        --sim_label $sim_label \
        --method $method_id \
        --sim_rds $sim_rds \
        --assignments_mtx $result_mtx \
        --out_prefix ${sim_label}__${method_id}
    """
}

process combine_assignments {
    label "split_combine"
    publishDir params.outdir, mode: 'copy'

    input:
    path sample_sheet
    path x

    output:
    path 'combined_results.Rds'

    script:
    """
    combine_assignments.R combined $sample_sheet $x
    """
}

process merge_confusion {
    label "split_combine"
    publishDir params.outdir, mode: 'copy'

    input:
    path confusion_list

    output:
    path "combined_confusion.csv"
    path 'combined_precision.png'
    path 'combined_recall.png'

    script:
    """
    combine_confusion.R combined $confusion_list
    """
}

// process and combine test stats
// FIXME: Reduce code duplication with process_assignments, combine_assignments?

process process_teststats {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method_id), path(sim_rds), path(result_mtx)

    output:
    path "${sim_label}__${method_id}_teststats.Rds"

    script:
    """
    process_teststats.R \
        --sim_label $sim_label \
        --method $method_id \
        --sim_rds $sim_rds \
        --stats_mtx $result_mtx \
        --out_rds ${sim_label}__${method_id}_teststats.Rds
    """
}

process combine_teststats {
    label "split_combine"
    publishDir params.outdir, mode: 'copy'

    input:
    path sample_sheet
    path x

    output:
    path 'combined_teststats.Rds'

    script:
    """
    combine_teststats.R combined_teststats.Rds $sample_sheet $x
    """
}
