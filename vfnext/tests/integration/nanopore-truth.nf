nextflow.enable.dsl = 2

include { NANOPORE } from '../../workflows/NANOPORE.nf'

workflow NANOPORE_TRUTH {
    take:
        readsCh
        ref

    main:
        NANOPORE(readsCh, ref)

    emit:
        filtered = NANOPORE.out.filteredVcfsCh
        consensus = NANOPORE.out.consensusCh
        summary = NANOPORE.out.summaryCh
}
