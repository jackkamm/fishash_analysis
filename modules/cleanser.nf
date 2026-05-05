process run_cleanser {
    label "cleanser"
    tag "seqtech=${seqtech},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_mtx)
    each seqtech

    output:
    tuple(val(sim_label),
          val(seqtech),
          path(sim_rds),
          path("posterior.mtx"))

    script:
    """
    cleanser -i $sim_mtx -o posterior.mtx --${seqtech} -p 8
    """
}

process convert_cleanser_out {
    label "data_conversion"

    input:
    tuple val(sim_label), val(seqtech), path(sim_rds), path(posterior_mtx)
    each cutoff

    output:
    tuple(val(sim_label),
          val("cleanser_${seqtech}_${cutoff}"),
          path(sim_rds),
          path("converted.mtx"))

    script:
    """
    convert_cleanser_out.R \
        --cleanser_posterior $posterior_mtx \
        --out_mtx converted.mtx \
        --cutoff $cutoff
    """
}
