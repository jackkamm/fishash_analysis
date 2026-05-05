include { run_fishash } from '../modules/fishash'

workflow FISHASH {
    take:
    sims_ch
    refit_list
    padj_list

    main:
    run_fishash(sims_ch, refit_list, padj_list)

    emit:
    out = run_fishash.out
}
