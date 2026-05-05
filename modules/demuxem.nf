process run_demuxem {
    label "demuxem"
    tag "minsignal=${min_signal},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)
    each min_signal

    output:
    tuple(val(sim_label),
          val("demuxem_signal${min_signal}"),
          path(sim_rds),
          path("out.mtx"))

    script:
    """
    run_demuxem.py $sim_h5ad out.mtx \
        --min_signal $min_signal --n_threads $task.cpus
    """
}
