include { get_gex_matrix_rds; run_sceptre_mixture } from '../modules/sceptre'

workflow SCEPTRE {
    take:
    sims_ch

    main:
    get_gex_matrix_rds(sims_ch)
    run_sceptre_mixture(get_gex_matrix_rds.out)

    emit:
    out = run_sceptre_mixture.out
}
