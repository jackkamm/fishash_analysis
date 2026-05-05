process convert_anndata {
    label 'data_conversion'

    input:
    tuple val(sim_label), path(sim_rds)

    output:
    tuple val(sim_label), path(sim_rds), path('anndata.h5ad')

    script:
    """
    convert_to_anndata.R \
        --in_rds $sim_rds \
        --out_h5ad anndata.h5ad
    """
}

process convert_mtx {
    label 'data_conversion'

    input:
    tuple val(sim_label), path(sim_rds)

    output:
    tuple val(sim_label), path(sim_rds), path('counts.mtx')

    script:
    """
    convert_to_mtx.R \
        --in_rds $sim_rds \
        --out_mtx counts.mtx
    """
}
