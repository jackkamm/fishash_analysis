include { convert_anndata }                                     from './modules/convert'
include { process_assignments; combine_assignments; merge_confusion; combine_teststats; process_teststats } from './modules/aggregate'
include { FISHASH }  from './workflows/fishash'
include { DEMUXEM }  from './workflows/demuxem'
include { CLEANSER } from './workflows/cleanser'
include { GEOMUX }   from './workflows/geomux'
include { SCEPTRE }  from './workflows/sceptre'

// HACK for rerunning fishash to get the continuous test statistics out
include { run_fishash as rerun_fishash } from './workflows/fishash'

workflow {
    samples_split = Channel.fromPath(params.sample_sheet).splitCsv(header:true)
    sims_ch = samples_split.map{ row -> tuple(row.sim_label, file(row.path))}

    precision_sweep = [.01, .05, .1, .2, .35, .5, .65, .8, .9, .95, .99]

    out_ch = channel.empty()
    teststats_ch = channel.empty()

    FISHASH(sims_ch, [10], precision_sweep)
    out_ch = out_ch.mix(FISHASH.out.out.map{
        simlab, methodlab, sim_rds, out_mtx, teststats_mtx ->
        tuple(simlab, methodlab, sim_rds, out_mtx)
    })

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
        out_ch = out_ch.mix(CLEANSER.out.assignments)
        teststats_ch = teststats_ch.mix(CLEANSER.out.stats)
    }

    process_assignments(out_ch)

    merge_confusion(process_assignments.out.confusion.collect())

    // rerun fishash at a single cutoff to get the continuous test statistics
    rerun_fishash(sims_ch, [10], [.05])
    teststats_ch = teststats_ch.mix(rerun_fishash.out.map{
        simlab, methodlab, sim_rds, out_mtx, teststats_mtx ->
        tuple(simlab, methodlab, sim_rds, teststats_mtx)
    })

    process_teststats(teststats_ch)

    combine_teststats(
        file(params.sample_sheet),
        process_teststats.out.collect()
    )

    // These can get large in this setting, and not really needed for plotting.
    //combine_assignments(
    //    file(params.sample_sheet),
    //    process_assignments.out.matrix.collect()
    //)
}
