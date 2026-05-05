process run_fishash {
    label "fishash"
    tag "refit=${refit_times},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds)
    each refit_times

    output:
    tuple(val(sim_label),
          val("fishash_refit${refit_times}"),
          path(sim_rds),
          path("out.mtx"))

    script:
    """
    run_fishash.R \
        --in_rds $sim_rds \
        --refit $refit_times \
        --out_mtx out.mtx
    """
}
