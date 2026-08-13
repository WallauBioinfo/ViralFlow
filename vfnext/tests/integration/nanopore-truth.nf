nextflow.enable.dsl = 2

include { NANOPORE } from '../../workflows/NANOPORE.nf'

workflow NANOPORE_TRUTH {
    take:
        reads_ch
        ref

    main:
        NANOPORE(reads_ch, ref)

    emit:
        filtered = NANOPORE.out.filtered_vcfs_ch
        consensus = NANOPORE.out.consensus_ch
        qc = NANOPORE.out.qc_ch
}
