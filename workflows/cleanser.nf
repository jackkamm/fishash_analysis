include { convert_mtx } from '../modules/convert'
include { run_cleanser; reformat_cleanser_mtx; convert_cleanser_out } from '../modules/cleanser'

workflow CLEANSER {
    take:
    sims_ch
    seqtech_list
    cutoff_list

    main:
    convert_mtx(sims_ch)
    run_cleanser(convert_mtx.out, seqtech_list)
    reformat_cleanser_mtx(run_cleanser.out)
    convert_cleanser_out(run_cleanser.out, cutoff_list)

    emit:
    assignments = convert_cleanser_out.out
    stats = reformat_cleanser_mtx.out
}
