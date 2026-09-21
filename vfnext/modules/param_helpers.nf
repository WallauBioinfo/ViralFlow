// Helpers for interpreting pipeline parameters.
//
// Separate from the workflows and from metadata_helpers.nf because both read
// the same values and must agree: NANOPORE decides whether to run the trimming
// step from exactly the value metadata_helpers uses to decide whether to record
// the tool that performs it. Two copies of the rule would let a run trim reads
// while its software_versions.tsv reported no trimming tool, or the reverse.

// nextflow.config declares trimLen as an Integer, but a value given on the
// command line arrives as a String, and the wrapper forwards params-file
// entries as command line arguments too. Comparing the two aborts the run
// before any task is submitted ("Cannot compare java.lang.String with value
// '30' and java.lang.Integer with value '0'"), so normalize once, here rather
// than in step0 validation: the integration tests call NANOPORE directly,
// without processInputs.
def normalizeTrimLen(value) {
    def raw = value == null ? '0' : value.toString().trim()
    if (!(raw ==~ /\d+/)) {
        error "--trimLen must be a non-negative integer, got '${value}'"
    }
    return raw as Integer
}
