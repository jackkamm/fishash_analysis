include { convert_anndata }                                     from './modules/convert'
include { process_assignments; combine_assignments; merge_confusion } from './modules/aggregate'
include { FISHASH }  from './workflows/fishash'
include { DEMUXEM }  from './workflows/demuxem'
include { CLEANSER } from './workflows/cleanser'
include { GEOMUX }   from './workflows/geomux'
include { SCEPTRE }  from './workflows/sceptre'

workflow {
    samples_split = Channel.fromPath(params.sample_sheet).splitCsv(header:true)
    sims_ch = samples_split.map{ row -> tuple(row.sim_label, file(row.path))}

    precision_sweep = [.01, .05, .1, .2, .35, .5, .65, .8, .9, .95, .99]

    out_ch = channel.empty()

    FISHASH(sims_ch, [10], precision_sweep)
    out_ch = out_ch.mix(FISHASH.out.out)

    if (!params.skipSceptre) {
        SCEPTRE(sims_ch, precision_sweep)
        out_ch = out_ch.mix(SCEPTRE.out.out)
    }

    convert_anndata(sims_ch)
    adata_ch = convert_anndata.out

    GEOMUX(adata_ch, [5], precision_sweep)
    out_ch = out_ch.mix(GEOMUX.out.out)

    DEMUXEM(adata_ch, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20])
    out_ch = out_ch.mix(DEMUXEM.out.out)

    if (!params.skipCleanser) {
        CLEANSER(sims_ch, ["cs", "dc"], precision_sweep)
        out_ch = out_ch.mix(CLEANSER.out.out)
    }

    process_assignments(out_ch)

    merge_confusion(process_assignments.out.confusion.collect())

    // These can get large in this setting, and not really needed for plotting.
    //combine_assignments(
    //    file(params.sample_sheet),
    //    process_assignments.out.matrix.collect()
    //)
}
