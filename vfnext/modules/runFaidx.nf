
process runFaidx {
    label "NP_basecontainer"
    tag "${ref}"

    input:
        path(ref)

    output:
        path("${ref}.fai"), emit: fai

    script:
    """
    set -euo pipefail

    samtools faidx ${ref}
    """
}

/*
Indexes the reference once per run, rather than inside every Clair3 task.

Deliberately separate from modules/genFaIdx.nf, which does the same thing for
ILLUMINA: configs/containers.config maps genFaIdx to the generate_consensus SIF,
which the nanopore container set does not build, so sharing it would fail under
-profile docker and for anyone who has only built the nanopore base container.
*/
