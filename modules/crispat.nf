process run_crispat_gauss {
    label "crispat_small"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_gauss"),
          path(sim_rds),
          path('crispat_gauss_assignments.csv'))

    script:
    """
    run_crispat_gauss.py $sim_h5ad crispat_gauss_
    rm -rf crispat_gauss_batch0/loss_plots
    """
}

process run_crispat_poisgauss {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_poisgauss"),
          path(sim_rds),
          path('crispat_poisgauss_assignments.csv'))

    script:
    """
    run_crispat_poisgauss.py $sim_h5ad crispat_poisgauss_
    rm -rf crispat_poisgauss_loss_plots
    """
}

process run_crispat_poisson {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_poisson"),
          path(sim_rds),
          path('crispat_poisson_assignments.csv'))

    script:
    """
    run_crispat_poisson.py $sim_h5ad crispat_poisson_ --cpus $task.cpus
    rm -rf crispat_poisson_loss_plots
    """
}

process run_crispat_negbinom {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_negbinom"),
          path(sim_rds),
          path('crispat_negbinom_assignments.csv'))

    script:
    """
    run_crispat_negbinom.py $sim_h5ad crispat_negbinom_ --cpus $task.cpus
    rm -rf crispat_negbinom_loss_plots
    """
}

process convert_crispat {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method), path(sim_rds), path(assignments_csv)

    output:
    tuple(val(sim_label),
          val(method),
          path(sim_rds),
          path("${method}__${sim_label}.mtx"))

    script:
    """
    convert_crispat.R \
        --sim_rds $sim_rds \
        --crispat_out $assignments_csv \
        --out_mtx ${method}__${sim_label}.mtx
    """
}
