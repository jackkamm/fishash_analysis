include { get_gex_matrix_rds; run_sceptre_mixture } from '../modules/sceptre'

workflow SCEPTRE {
    take:
    sims_ch
    prob_list

    main:
    get_gex_matrix_rds(sims_ch)
    run_sceptre_mixture(get_gex_matrix_rds.out, prob_list)

    emit:
    out = run_sceptre_mixture.out
}
