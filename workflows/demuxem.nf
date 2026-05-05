include { run_demuxem } from '../modules/demuxem'

workflow DEMUXEM {
    take:
    adata_ch
    min_signal_list

    main:
    run_demuxem(adata_ch, min_signal_list)

    emit:
    out = run_demuxem.out
}
