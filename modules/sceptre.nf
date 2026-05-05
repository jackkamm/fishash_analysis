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
    tag "simlab=${sim_label},prob=${probability_threshold}"

    input:
    tuple val(sim_label), path(sim_rds), path(gex_rds)
    each probability_threshold

    output:
    tuple(val(sim_label),
          val("sceptre_mixture_prob${probability_threshold}"),
          path(sim_rds),
          path("sceptre.mtx"))

    script:
    """
    run_sceptre_mixture.R \
        --in_rds $sim_rds \
        --gex_rds $gex_rds \
        --out_mtx sceptre.mtx \
        --probability_threshold $probability_threshold \
        --cpus $task.cpus
    """
}
