include { run_geomux; convert_geomux } from '../modules/geomux'

workflow GEOMUX {
    take:
    adata_ch
    min_umi_list

    main:
    run_geomux(adata_ch, min_umi_list)
    convert_geomux(run_geomux.out)

    emit:
    out = convert_geomux.out
}
