process run_fishash {
    label "fishash"
    tag "refit=${refit_times},simlab=${sim_label},padj=${padj_cutoff}"

    input:
    tuple val(sim_label), path(sim_rds)
    each refit_times
    each padj_cutoff

    output:
    tuple(val(sim_label),
          val("fishash_refit${refit_times}_padj${padj_cutoff}"),
          path(sim_rds),
          path("out.mtx"))

    script:
    """
    run_fishash.R \
        --in_rds $sim_rds \
        --refit $refit_times \
        --padj_cutoff $padj_cutoff \
        --out_mtx out.mtx
    """
}
