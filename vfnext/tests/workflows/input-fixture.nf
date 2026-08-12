#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { processInputs } from '../../workflows/step0-input-handling.nf'

workflow {
    processInputs()
    processInputs.out.reads_ch.view { meta, reads ->
        "PREPARED ${meta.id} paired=${meta.is_paired_end} files=${reads*.name.join(',')} batch=${meta.batch ?: ''}"
    }
}
