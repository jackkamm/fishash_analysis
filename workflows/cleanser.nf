include { convert_mtx } from '../modules/convert'
include { run_cleanser; convert_cleanser_out } from '../modules/cleanser'

workflow CLEANSER {
    take:
    sims_ch
    seqtech_list
    cutoff_list

    main:
    convert_mtx(sims_ch)
    run_cleanser(convert_mtx.out, seqtech_list)
    convert_cleanser_out(run_cleanser.out, cutoff_list)

    emit:
    out = convert_cleanser_out.out
}
