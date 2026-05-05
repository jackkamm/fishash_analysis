include { run_crispat_gauss;
          run_crispat_poisgauss;
          run_crispat_poisson;
          run_crispat_negbinom;
          convert_crispat } from '../modules/crispat'

// NOTE: A difficulty with crispat is that sometimes the dask
// workers give timeout errors. It helps to run the pipeline with
// retries. Possibly should modify the timeout params within the
// crispat implementations to avoid this

workflow CRISPAT {
    take:
    adata_ch

    main:
    run_crispat_gauss(adata_ch)
    crispat_out = run_crispat_gauss.out

    run_crispat_poisgauss(adata_ch)
    crispat_out = crispat_out.mix(run_crispat_poisgauss.out)

    if (!params.skipCrispatPoisson) {
        run_crispat_poisson(adata_ch)
        crispat_out = crispat_out.mix(run_crispat_poisson.out)
    }

    if (!params.skipCrispatNegBinom) {
        run_crispat_negbinom(adata_ch)
        crispat_out = crispat_out.mix(run_crispat_negbinom.out)
    }

    convert_crispat(crispat_out)

    emit:
    out = convert_crispat.out
}
