process run_geomux {
    label "fishash"
    tag "minumi=${min_umi},simlab=${sim_label},fdr=${fdr_threshold}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)
    each min_umi
    each fdr_threshold

    output:
    tuple(val(sim_label),
          val("geomux_minumi${min_umi}_fdr${fdr_threshold}"),
          path(sim_rds),
          path("geomux_out.tsv"))

    script:
    """
    geomux \
        $sim_h5ad \
        geomux_out.tsv \
        --fdr-threshold $fdr_threshold \
        --min-umi-cells ${min_umi} \
        --min-umi-guides ${min_umi}
    """
}

process convert_geomux {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method), path(sim_rds), path(geomux_tsv)

    output:
    tuple(val(sim_label),
          val(method),
          path(sim_rds),
          path("geomux_converted.mtx"))

    script:
    """
    convert_geomux_tsv.R \
        --geomux_tsv $geomux_tsv \
        --orig_rds $sim_rds \
        --out_mtx geomux_converted.mtx
    """
}
