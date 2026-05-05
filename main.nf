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

process run_fishash {
    label "fishash"
    tag "refit=${refit_times},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds)
    each refit_times

    output:
    tuple(val(sim_label),
          val("fishash_refit${refit_times}"),
          path(sim_rds),
          path("out.mtx"))

    script:
    """
    run_fishash.R \
        --in_rds $sim_rds \
        --refit $refit_times \
        --out_mtx out.mtx
    """
}

process run_demuxem {
    label "demuxem"
    tag "minsignal=${min_signal},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)
    each min_signal

    output:
    tuple(val(sim_label),
          val("demuxem_signal${min_signal}"),
          path(sim_rds),
          path("out.mtx"))

    script:
    """
    run_demuxem.py $sim_h5ad out.mtx \
        --min_signal $min_signal --n_threads $task.cpus
    """
}

process run_cleanser {
    label "cleanser"
    tag "seqtech=${seqtech},simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_mtx)
    each seqtech

    output:
    tuple(val(sim_label),
          val(seqtech),
          path(sim_rds),
          path("posterior.mtx"))

    script:
    """
    cleanser -i $sim_mtx -o posterior.mtx --${seqtech} -p 8
    """
}

process convert_cleanser_out {
    label "data_conversion"

    input:
    tuple val(sim_label), val(seqtech), path(sim_rds), path(posterior_mtx)
    each cutoff

    output:
    tuple(val(sim_label),
          val("cleanser_${seqtech}_${cutoff}"),
          path(sim_rds),
          path("converted.mtx"))

    script:
    """
    convert_cleanser_out.R \
        --cleanser_posterior $posterior_mtx \
        --out_mtx converted.mtx \
        --cutoff $cutoff
    """
}

process run_geomux {
    label "fishash"
    tag "minumi=${min_umi},simlab=${sim_label}"
    
    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)
    each min_umi

    output:
    tuple(val(sim_label),
          val("geomux_minumi${min_umi}"),
          path(sim_rds),
          path("geomux_out.tsv"))

    script:
    """
    geomux \
        $sim_h5ad \
        geomux_out.tsv \
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

process convert_crispat {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method), path(sim_rds), path(assignments_csv)

    output:
    tuple(val(sim_label),
          val(method),
          path(sim_rds),
          path("${method}__${sim_label}.mtx"))

    script:
    """
    convert_crispat.R \
        --sim_rds $sim_rds \
        --crispat_out $assignments_csv \
        --out_mtx ${method}__${sim_label}.mtx
    """
}

process run_crispat_gauss {
    label "crispat_small"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_gauss"),
          path(sim_rds),
          path('crispat_gauss_assignments.csv'))

    script:
    """
    run_crispat_gauss.py $sim_h5ad crispat_gauss_
    """
}

process run_crispat_poisgauss {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_poisgauss"),
          path(sim_rds),
          path('crispat_poisgauss_assignments.csv'))

    script:
    """
    run_crispat_poisgauss.py $sim_h5ad crispat_poisgauss_
    """
}

process run_crispat_poisson {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_poisson"),
          path(sim_rds),
          path('crispat_poisson_assignments.csv'))

    script:
    """
    run_crispat_poisson.py $sim_h5ad crispat_poisson_ --cpus $task.cpus
    """
}

process run_crispat_negbinom {
    label "crispat_big"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(sim_h5ad)

    output:
    tuple(val(sim_label),
          val("crispat_negbinom"),
          path(sim_rds),
          path('crispat_negbinom_assignments.csv'))

    script:
    """
    run_crispat_negbinom.py $sim_h5ad crispat_negbinom_ --cpus $task.cpus
    """
}

process get_gex_matrix_rds {
    label "sceptre"

    input:
    tuple val(sim_label), path(sim_rds)

    output:
    tuple(val(sim_label),
          path(sim_rds),
          path("gex.Rds"))

    script:
    """
    get_gex_matrix.R \
        --in_rds $sim_rds \
        --out_rds gex.Rds
    """
}

process run_sceptre_mixture {
    label "sceptre"
    tag "simlab=${sim_label}"

    input:
    tuple val(sim_label), path(sim_rds), path(gex_rds)

    output:
    tuple(val(sim_label),
          val("sceptre_mixture"),
          path(sim_rds),
          path("sceptre.mtx"))

    script:
    """
    run_sceptre_mixture.R \
        --in_rds $sim_rds \
        --gex_rds $gex_rds \
        --out_mtx sceptre.mtx \
        --cpus $task.cpus
    """
}

// process individual assignment results to prepare for merging them
// computes confusion matrices, and adds metadata to assignment matrices
process process_assignments {
    label "data_conversion"

    input:
    tuple val(sim_label), val(method_id), path(sim_rds), path(result_mtx)

    output:
    path "${sim_label}__${method_id}_confusion.csv", emit: confusion
    path "${sim_label}__${method_id}_matrixWithMeta.Rds", emit: matrix

    script:
    """
    process_assignments.R \
        --sim_label $sim_label \
        --method $method_id \
        --sim_rds $sim_rds \
        --assignments_mtx $result_mtx \
        --out_prefix ${sim_label}__${method_id}
    """
}

process combine_assignments {
    label "split_combine"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path sample_sheet
    path x

    output:
    path 'combined_results.Rds'

    script:
    """
    combine_assignments.R combined $sample_sheet $x
    """
}

process merge_confusion {
    label "split_combine"
    publishDir params.outdir, mode: 'copy'

    input:
    path confusion_list

    output:
    path "combined_confusion.csv"
    path 'combined_precision.png'
    path 'combined_recall.png'

    script:
    """
    combine_confusion.R combined $confusion_list
    """
}

workflow {
    samples_split = Channel.fromPath(params.sample_sheet).splitCsv(header:true)
    sims_ch = samples_split.map{ row -> tuple(row.sim_label, file(row.path))}

    out_ch = channel.empty()

    refit_times = [0, 10]
    run_fishash(sims_ch, refit_times)
    out_ch = out_ch.mix(run_fishash.out)

    if (!params.skipSceptre) {
        get_gex_matrix_rds(sims_ch)
        run_sceptre_mixture(get_gex_matrix_rds.out)
        out_ch = out_ch.mix(run_sceptre_mixture.out)
    }

    convert_anndata(sims_ch)
    adata_ch = convert_anndata.out

    run_geomux(adata_ch, [1, 5])
    convert_geomux(run_geomux.out)
    out_ch = out_ch.mix(convert_geomux.out)

    min_signal = [2, 10]
    run_demuxem(adata_ch, min_signal)
    out_ch = out_ch.mix(run_demuxem.out)

    if (!params.skipCleanser) {
        convert_mtx(sims_ch)
        mtx_ch = convert_mtx.out

        run_cleanser(mtx_ch, ["cs", "dc"])
        convert_cleanser_out(run_cleanser.out, [0.5, 0.95])
        out_ch = out_ch.mix(convert_cleanser_out.out)
    }

    run_crispat_gauss(adata_ch)
    crispat_out = run_crispat_gauss.out

    run_crispat_poisgauss(adata_ch)
    crispat_out = crispat_out.mix(run_crispat_poisgauss.out)

    // NOTE: A difficulty with crispat is that sometimes the dask
    // workers give timeout errors. It helps to run the pipeline with
    // retries. Possibly should modify the timeout params within the
    // crispat implementations to avoid this

    if (!params.skipCrispatPoisson) {
        run_crispat_poisson(adata_ch)
        crispat_out = crispat_out.mix(run_crispat_poisson.out)
    }

    if (!params.skipCrispatNegBinom) {
        run_crispat_negbinom(adata_ch)
        crispat_out = crispat_out.mix(run_crispat_negbinom.out)
    }

    convert_crispat(crispat_out)

    out_ch = out_ch.mix(convert_crispat.out)

    process_assignments(out_ch)

    merge_confusion(process_assignments.out.confusion.collect())

    combine_assignments(
        file(params.sample_sheet),
        process_assignments.out.matrix.collect()
    )
}
