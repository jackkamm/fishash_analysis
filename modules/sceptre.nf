process get_gex_matrix_rds {
    label "sceptre"

    input:
    tuple val(sim_label), path(sim_rds)

    output:
    tuple(val(sim_label),
          path(sim_rds),
          path("gex.Rds"))

    script:
    """
    get_gex_matrix.R \
        --in_rds $sim_rds \
        --out_rds gex.Rds
    """
}

process run_sceptre_mixture {
    label "sceptre"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(gex_rds)

    output:
    tuple(val(sim_label),
          val("sceptre_mixture"),
          path(sim_rds),
          path("sceptre.mtx"))

    script:
    """
    run_sceptre_mixture.R \
        --in_rds $sim_rds \
        --gex_rds $gex_rds \
        --out_mtx sceptre.mtx \
        --cpus $task.cpus
    """
}
