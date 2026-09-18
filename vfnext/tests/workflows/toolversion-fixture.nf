nextflow.enable.dsl = 2

include { capture_tool_version } from '../../modules/metadata.nf'

// Drives capture_tool_version with `echo` standing in for each tool, so the
// version parsing can be pinned against the formats the pipeline's tools really
// print without needing those tools present. The strings below were captured
// from live runs (clair3, porechop_abi, minimap2, samtools, bcftools) or are the
// documented output of the ILLUMINA tools.
//
// A greedy sed expression used to anchor on the last digit group in the line, so
// everything with three components was recorded wrong - "Clair3 v1.2.0" as
// "2.0", "fastp 0.23.2" as "3.2". Only two-component versions survived, and
// bcftools was the only case any test covered.
workflow TOOLVERSION_FIXTURE {
    main:
        def container = params.base_container.toString()

        channel.of(
            tuple("NANOPORE", "clair3", 'echo "Clair3 v1.2.0"', container),
            tuple("NANOPORE", "porechop_abi", 'echo "0.5.1"', container),
            tuple("NANOPORE", "minimap2", 'echo "2.28-r1209"', container),
            tuple("NANOPORE", "bcftools", 'echo "bcftools 1.21"', container),
            tuple("ILLUMINA", "fastp", 'echo "fastp 0.23.2"', container),
            tuple("ILLUMINA", "ivar", 'echo "iVar version 1.3.1"', container),
            tuple("ILLUMINA", "mafft", 'echo "v7.505 (2022/Apr/10)"', container),
            tuple("ILLUMINA", "nextclade", 'echo "nextclade 3.18.0"', container),
            tuple("ILLUMINA", "unversioned", 'echo "no version here"', container)
        )
        .set { specs_ch }

        capture_tool_version(specs_ch)

    emit:
        capture_tool_version.out
}
